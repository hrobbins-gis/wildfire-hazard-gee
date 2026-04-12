// -----------------------------
// WILDFIRE HAZARD MODEL
// -----------------------------

// You can manually draw or put in your vertices
// make sure your polygon is named "aoi"
// var aoi =
// START AND END ATECEDENT TIME PERIOD OF FIRE
var START = ee.Date('2023-10-01');
var END = ee.Date('2024-09-01');

// CHANGE VARIABLE WEIGHTS
var weights = {
  precip: 0.30,   // sum of precipitation
  vs: 0.10,  // mean wind speed velocity
  vpd: 0.15,  // mean vapor pressure deficit
  ndvi: 0.10,   // min NDVI
  slope: 0.35   // slope
};

// EXPORT FOLDER SET UP - CHANGE TO YOUR PREFERRED FOLDER IN GOOGLE DRIVE
var exportFolder = 'YourFolderName';

// -------------------------------------------
// HELPER FUNCTIONS - DO NOT CHANGE + COLLAPSE
// -------------------------------------------

// Cloud mask function
function maskSrClouds(image) {
var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111',2)).eq(0);
var saturationMask = image.select('QA_RADSAT').eq(0);
return image.updateMask(qaMask).updateMask(saturationMask);
}

// Robust extractor of a 6-element affine transform from any projection
// Needed for reprojecting before export and accurate slope calculation
// Adapted by ChatGPT
function getTransformListFromProjection(proj) {
  // Use getInfo() to pull the transform metadata to the client
  var info = proj.transform().getInfo(); // safe: tiny object

  // Case 1: info is an object with a 'matrix' key (3x3)
  if (info && typeof info === 'object' && info.matrix) {
    var m = info.matrix;
    return [
      m[0][0], // a
      m[0][1], // b
      m[0][2], // c
      m[1][0], // d
      m[1][1], // e
      m[1][2]  // f
    ];
  }

  // Case 2: info is an object with scale/translate (and optional shear)
  if (info && typeof info === 'object' && info.scale && info.translate) {
    var shear = info.shear || [0, 0];
    return [
      info.scale[0],
      shear[0],
      info.translate[0],
      shear[1],
      info.scale[1],
      info.translate[1]
    ];
  }

  // Case 3: info is a string (PARAM_MT). Parse only the elt_* parameters.
  if (typeof info === 'string') {
    var txt = info;
    // Regex to find PARAMETER["elt_i_j", value]
    var re = /PARAMETER\["elt_(\d)_(\d)",\s*(-?\d+\.?\d*)\]/g;
    var vals = {
      "0_0": 0, "0_1": 0, "0_2": 0,
      "1_0": 0, "1_1": 0, "1_2": 0
    };
    var m;
    while ((m = re.exec(txt)) !== null) {
      var key = m[1] + "_" + m[2];     // e.g. "0_2"
      var v = parseFloat(m[3]);
      vals[key] = v;
    }
    // If we didn't find any elt_ parameters, try a looser numeric fallback:
    var foundAny = Object.keys(vals).some(function(k){ return vals[k] !== 0; });
    if (!foundAny) {
      // try to grab numbers ignoring metadata labels (last resort)
      var allNums = txt.match(/-?\d+\.?\d*/g);
      if (allNums && allNums.length >= 6) {
        // assume the last 6 numbers are the affine elements (risky fallback)
        var L = allNums.length;
        return [
          parseFloat(allNums[L-6]),
          parseFloat(allNums[L-5]),
          parseFloat(allNums[L-4]),
          parseFloat(allNums[L-3]),
          parseFloat(allNums[L-2]),
          parseFloat(allNums[L-1])
        ];
      }
    }
    return [
      vals["0_0"], vals["0_1"], vals["0_2"],
      vals["1_0"], vals["1_1"], vals["1_2"]
    ];
  }

  // If we get here, we couldn't parse it
  throw new Error('Could not parse transform from projection. getInfo() returned: ' + JSON.stringify(info));
}

// Calculate NDVI function with scaling surface reflectance
function addNDVI(image) {
  // Scale reflectance
  var sr = image.multiply(0.0000275).add(-0.2);
  // NDVI calculation: (NIR - RED) / (NIR + RED)
  var ndvi = sr.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  return sr.addBands(ndvi).copyProperties(image, image.propertyNames());
}

// Calculate VPD function
var addVPD = function(img) {
  var tmin = img.select('tmin');
  var tmax = img.select('tmax');
  var vp = img.select('vp').divide(1000); // already in kPa

  // Saturation vapor pressure from tmin
  var es_tmin = tmin.expression(
    '0.6108 * exp((17.27 * T) / (T + 237.3))', {
      'T': tmin
    });

  // Saturation vapor pressure from tmax
  var es_tmax = tmax.expression(
    '0.6108 * exp((17.27 * T) / (T + 237.3))', {
      'T': tmax
    });

  // Mean saturation vapor pressure
  var es_mean = es_tmin.add(es_tmax).divide(2);

  // VPD = es_mean - actual vapor pressure
  var vpd = es_mean.subtract(vp)
    .max(0)
    .rename('vpd');

  return img.addBands(vpd);
};

// MIN-MAX NORMALIZE (per AOI) ----------
function minMaxNormalize(img, bandName, region) {
  // compute min/max within region
  var stats = img.select(bandName).reduceRegion({
    reducer: ee.Reducer.percentile([0,100]),
    geometry: region,
    scale: 30,
    bestEffort: true,
    maxPixels: 1e13
  });
  var minVal = ee.Number(stats.get(bandName + '_p0'));
  var maxVal = ee.Number(stats.get(bandName + '_p100'));
  // avoid divide by zero
  var denom = maxVal.subtract(minVal).max(1e-6);
  return img.select(bandName).subtract(minVal).divide(denom).rename(bandName + '_norm');
}

// Inspect image function
// Prints CRS, scale and minimum/maximum values
// Use for checking intermediate rasters
function inspectImage(img, name) {

  print('--- ' + name + ' ---');

  // Projection (safe)
  print(name + ' CRS:', img.projection().crs());
  print(name + ' nominal scale:', img.projection().nominalScale());

  // Min/max sampled safely using a small, reduced region
  // Uses bestEffort + a coarser scale to avoid heavy operations
  var stats = img.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: img.geometry().bounds(),  // bounded region
    scale: 200,                         // coarser for speed
    bestEffort: true,
    maxPixels: 1e9
  });

  print(name + ' min/max:', stats);

  print('----------------------------');
}

// ---- AOI SIZE CHECK ----
var aoiAreaHa = aoi.area(1,null).divide(10000);  // m² → hectares
var aoiAreaKm2 = aoiAreaHa.divide(100);    // hectares → km²

print('AOI area (km²):', aoiAreaKm2);

// Tested threshold - at 30 m resolution, any larger aoi will likely time out
var maxRecommendedKm2 = 2000;

// Large AOI warning message:
if (aoiAreaKm2.gt(maxRecommendedKm2).getInfo()) {
  print('⚠️ WARNING: AOI is large (' + aoiAreaKm2.getInfo().toFixed(1) +
        ' km²). \nSome layers and calculations may time out. \nConsider using a smaller analysis polygon.');
}

// Pull in World Land Cover
var wc = ee.ImageCollection("ESA/WorldCover/v100").first();

// Create burn mask
var burnable = wc.eq(10)  // forest
  .or(wc.eq(20))         // shrubland
  .or(wc.eq(30))          // grassland
  .or(wc.eq(40))         // crops
  .clip(aoi);

// Pull in Landsat 8 and 9 images
// Change cloud cover filter if desired
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');
var landsat = l8.merge(l9)
  .filterBounds(aoi)
  .filterDate(START, END)
  .filter(ee.Filter.lt('CLOUD_COVER', 15)); // optional

// Apply the cloud mask to the collection
var landsatFiltMasked = landsat.map(maskSrClouds);

// Variables for reference, projection and transformation used for exports and slope calc.
var refBand = landsat.first().select('SR_B4');
var refProj = refBand.projection();
var refCRS = refProj.crs().getInfo();
var refTran = getTransformListFromProjection(refProj);

// Uncomment for number of Landsat images in the console
// print("Cloud Masked Landsat Collection", landsatFiltMasked);

// Calculate NDVI over the image collection
var withNDVI = landsatFiltMasked.map(addNDVI).map(function(img) {return img.clip(aoi);});

// Compute the min NDVI composite (creates single image for NDVI)
// setDefaultProjection here converts to the Landsat projection, otherwise
// would default to EPSG 4326
var minNDVI = withNDVI.select('NDVI').min().setDefaultProjection(refProj).rename('ndvi').clip(aoi);

// Pull in DEM
var dem = ee.Image('USGS/SRTMGL1_003')
    .clip(aoi)
    .resample('bilinear')
    // reproject here for accurate slope units
    .reproject({
      crs: refCRS,
      crsTransform: refTran
    });

// Create slope image
var slope = ee.Terrain.slope(dem).rename('slope').clip(aoi);

// Pull in gridMET data for wind speed (m/s)
var gridmet = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')
  .filterBounds(aoi)
  .filterDate(START, END);

// Create wind speed image, calculate mean (average)
var wind = gridmet.select('vs').mean().clip(aoi).rename('wind_mean');

// Pull in Daymet data for precipitation (mm) and vapor pressure (kPa)
var daymet = ee.ImageCollection('NASA/ORNL/DAYMET_V4')
  .filterDate(START, END)
  .filterBounds(aoi);

// Create SUM precipitation image
var precip = daymet.select('prcp').sum().clip(aoi).rename('precip_total');

// Add daily VPD values to each Daymet image
var daymetWithVPD = daymet.map(addVPD);

// Create mean VPD image
var vpd_mean = daymetWithVPD.select('vpd').mean().rename('vpd_mean').clip(aoi);

// Apply the burn mask
var ndvi_masked   = minNDVI.updateMask(burnable);
var slope_masked  = slope.updateMask(burnable);
var wind_masked   = wind.updateMask(burnable);
var precip_masked = precip.updateMask(burnable);
var vpd_masked    = vpd_mean.updateMask(burnable);

// Normalize to common scale 0 - 1
var norm_ndvi = minMaxNormalize(ndvi_masked, 'ndvi', aoi).clamp(0, 1)
var norm_slope = minMaxNormalize(slope_masked, 'slope', aoi).clamp(0, 1)
var norm_wind = minMaxNormalize(wind_masked, 'wind_mean', aoi).clamp(0, 1)
var norm_precip = minMaxNormalize(precip_masked, 'precip_total', aoi).clamp(0, 1)
var norm_vpd = minMaxNormalize(vpd_masked, 'vpd_mean', aoi).clamp(0, 1)

// Uncomment to check original/normalized values in the console
// NDVI
// inspectImage(ndvi_masked, 'NDVI Unnormalized');
// inspectImage(norm_ndvi, 'NDVI Normalized');

// Slope
// inspectImage(slope_masked, 'Slope Unnormalized');
// inspectImage(norm_slope, 'Slope Normalized');

// Wind Speed
// inspectImage(wind_masked, 'Wind Unnormalized');
// inspectImage(norm_wind, 'VS Normalized');

// Precipitation
// inspectImage(precip_masked, 'Precip Unnormalized');
// inspectImage(norm_precip, 'Precip Normalized');

// Vapor Pressure Deficit
// inspectImage(vpd_masked, 'VPD Unnormalized');
// inspectImage(norm_vpd, 'VPD Normalized');

// Convert to hazard variables
var risk_wind = norm_wind.rename('wind_risk');
var risk_vpd = norm_vpd.rename('vpd_risk');
var risk_slope = norm_slope.rename('slope_risk');

// Invert precipitation/ndvi
var risk_precip = ee.Image(1).subtract(norm_precip)
  .rename('precip_risk');

var risk_ndvi = ee.Image(1).subtract(norm_ndvi)
  .rename('ndvi_risk');

// ---------- RISK STACK ----------
var riskStack = ee.Image.cat([
  risk_precip,
  risk_wind,
  risk_vpd,
  risk_ndvi,
  risk_slope
]);

// ---------- WEIGHTED SUM ----------
var hazard = riskStack.expression(
  'wp*p + ww*w + wv*v + wn*n + ws*s', {
    wp: weights.precip, p: risk_precip,
    ww: weights.vs,     w: risk_wind,
    wv: weights.vpd,    v: risk_vpd,
    wn: weights.ndvi,   n: risk_ndvi,
    ws: weights.slope,  s: risk_slope
  }
).rename('hazard_index')
.clip(aoi);

var norm_hazard = minMaxNormalize(hazard, 'hazard_index', aoi).clamp(0, 1)

// Uncomment to check hazard values
// inspectImage(hazard, 'Original Hazard');
// inspectImage(norm_hazard, 'Hazard Normalized');

// ---------- ADD TO MAP ----------
Map.centerObject(aoi, 10);
// Add layers to the map for display
// 3 - 4 is advised to avoid time outs
// Min/Max will need to be scaled dependent on the selected variables,
// refer to citation docs for ranges or stat print outs using InspectImage function

// -- PRECIPITATION --
// Map.addLayer(precip, {min: 456, max: 1024, palette: ['brown','orange', 'yellow', 'green']}, 'OG Precip');
// Map.addLayer(precip_masked, {min: 456, max: 1024, palette: ['brown','orange', 'yellow', 'green']}, 'Masked Precip');
// Map.addLayer(norm_precip, {min: 0, max: 1, palette: ['brown','orange', 'yellow', 'green']}, 'Normalized Precip');
// Map.addLayer(risk_precip, {min: 0, max: 1, palette: ['green', 'yellow', 'red']}, 'Risk Precip');

// -- WIND SPEED (V/S)
// Map.addLayer(wind, {min: 3.20, max: 3.58, palette: ['green', 'yellow', 'red']}, 'OG Wind')
// Map.addLayer(wind_masked, {min: 3.20, max: 3.58, palette: ['green', 'yellow', 'red']}, 'Masked Wind')
// Map.addLayer(norm_wind, {min: 0, max: 1, palette: ['green', 'yellow', 'red']}, 'Norm Wind')
// Map.addLayer(risk_wind, {min: 0, max: 1, palette: ['green', 'yellow', 'red']}, 'Risk Wind')


// -- VAPOR PRESSURE DEFICIT (KPA)
// Map.addLayer(vpd_mean, {min: 0.62, max: 1.58, palette: ['green', 'yellow', 'red']}, 'VPD masked');
// Map.addLayer(vpd_masked, {min: 0.62, max: 1.58, palette: ['green', 'yellow', 'red']}, 'VPD masked');
// Map.addLayer(norm_vpd, {min: 0, max: 1, palette: ['green', 'yellow', 'red']}, 'Normalized VPD');
// Map.addLayer(risk_vpd, {min: 0, max: 1, palette: ['green', 'yellow', 'red']}, 'Risk VPD');

// -- SLOPE
// Map.addLayer(slope, {min: 0, max: 44, palette: ['green', 'yellow', 'red']}, 'OG Slope');
// Map.addLayer(slope_masked, {min: 0, max: 44, palette: ['green', 'yellow', 'red']}, 'Masked Slope');
// Map.addLayer(norm_slope, {min: 0, max: 1, palette: ['green', 'yellow', 'red']}, 'Normalized Slope');
// Map.addLayer(risk_slope, {min: 0, max: 0.66, palette: ['green', 'yellow', 'red']}, 'Risk Slope');


// -- NDVI
// Map.addLayer(minNDVI, {min: -0.99, max: 0.798, palette: ['brown','orange', 'yellow', 'green']}, 'Min NDVI');
// Map.addLayer(ndvi_masked, {min: -0.68, max: 0.798, palette: ['brown','orange', 'yellow', 'green']}, 'Masked NDVI');
// Map.addLayer(norm_ndvi, {min: 0, max: 0.95, palette: ['green', 'yellow', 'red']}, 'Normalized NDVI');
// Map.addLayer(risk_ndvi, {min: 0, max: 0.83, palette: ['green', 'yellow', 'red']}, 'Risk NDVI');


// Map.addLayer(burnable, null, 'Burn Mask');
// Map.addLayer(landcover, null, 'Landcover');

// Map.addLayer(norm_vpd, {min: 0, max: 1, palette: ['green', 'yellow', 'red']}, 'Daymet Normalized VPD');
// Map.addLayer(risk_vpd_gm, {min: 0, max: 1, palette: ['green', 'yellow', 'red']}, 'Gridmet Normalized VPD');
// Map.addLayer(hazard, {min:0, max:1, palette: ['green', 'yellow', 'red']}, 'Hazard Index');
Map.addLayer(norm_hazard, {min:0, max:1, palette: ['green', 'yellow', 'red']}, 'Normalized Hazard Index');

// -----------------------------
// ---------- EXPORTS ----------
// -----------------------------
// Uncomment for exports to designated Google Drive folder -- line 28

// -------- BURN MASK ----------
// var burnable_mask = burnable.selfMask();
// Export.image.toDrive({
//   image: burnable_mask,
//   description: 'Burnable_Mask',
//   folder: exportFolder,
//   scale: 30,
//   region: aoi,
//   crs: refCRS,
//   maxPixels: 1e13
// });

//  ------ PRECIPITATION ------
// Export.image.toDrive({
//   image: precip,
//   description: 'precip_sum',
//   folder: exportFolder,
//   fileNamePrefix: 'precip_sum',
//   region: aoi,
//   scale: 30,
//   crs: refCRS,
//   maxPixels: 1e13
// });

//  ------- WIND SPEED -------
// Export.image.toDrive({
//   image: wind,
//   description: 'wind_speed',
//   folder: exportFolder,
//   fileNamePrefix: 'wind_speed',
//   region: aoi,
//   scale: 30,
//   crs: refCRS,
//   maxPixels: 1e13
// });

//  ---------- VPD ----------
// Export.image.toDrive({
//   image: vpd_mean,
//   description: 'VPD_Mean',
//   folder: exportFolder,
//   fileNamePrefix: 'Mean_VPD',
//   region: aoi,
//   scale: 30,
//   crs: refCRS,
//   maxPixels: 1e13
// });

//  ---------- NDVI ----------
// Export.image.toDrive({
//   image: minNDVI,
//   description: 'min_ndvi',
//   folder: exportFolder,
//   fileNamePrefix: 'Min_NDVI',
//   region: aoi,
//   scale: 30,
//   crs: refCRS,
//   maxPixels: 1e13
// }


//  ---------- SLOPE ----------
// Export.image.toDrive({
//   image: slope,
//   description: 'Slope',
//   folder: exportFolder,
//   fileNamePrefix: 'Slope',
//   region: aoi,
//   scale: 30,
//   crs: refCRS,
//   maxPixels: 1e13
// });

//  ---------- FIRE HAZARD ----------
// Export.image.toDrive({
//   image: hazard,
//   description: 'Wildfire_Hazard',
//   folder: exportFolder,
//   fileNamePrefix: 'wildfire_hazard',
//   region: aoi,
//   scale: 30,
//   crs: refCRS,
//   maxPixels: 1e13
// });


//  ------- NORMALIZED FIRE HAZARD -------
// Export.image.toDrive({
//   image: norm_hazard,
//   description: 'Norm Wildfire_Hazard',
//   folder: exportFolder,
//   fileNamePrefix: 'norm_wildfire_hazard',
//   region: aoi,
//   scale: 30,
//   crs: refCRS,
//   maxPixels: 1e13
// });

// PLOTS FOR THE CONSOLE ONLY
var hazHistogram = ui.Chart.image.histogram({
  image: hazard,
  region: aoi,
  scale: 30,
  maxPixels: 1e9
})
.setOptions({
  title: 'WHI Distribution',
  hAxis: {title: 'Hazard Value'},
  vAxis: {title: 'Pixel Count'},
  colors: ['#4caf50']
});
print(hazHistogram);