// Make sure your polygon is named 'aoi'

var exportFolder = '';

var pre_fireSTART = ee.Date('2024-07-01');
var pre_fireEND = ee.Date('2024-09-01');

var post_fireSTART = ee.Date('2024-09-04');
var post_fireEND = ee.Date('2025-01-31');

// Define the cloud mask function.
function maskSrClouds(image) {
var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111',2)).eq(0);
var saturationMask = image.select('QA_RADSAT').eq(0);
return image.updateMask(qaMask).updateMask(saturationMask);
}

// Function to calculate NBR
function addNBR(image) {

  // Apply reflectance scaling (Collection 2 L2)
  var sr = image.select(['SR_B.*'])
                .multiply(0.0000275)
                .add(-0.2);

  // Calculate NBR using NIR (B5) and SWIR2 (B7)
  var nbr = sr.normalizedDifference(['SR_B5', 'SR_B7'])
              .rename('NBR');

  return image.addBands(nbr);
}

// Pull in Landsat 8 and 9 images, select time frame and cloud cover threshold
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2');
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2');



var pre_fire = l8.merge(l9)
  .filterBounds(aoi)
  .filterDate(pre_fireSTART,pre_fireEND)
  .map(maskSrClouds)
  .map(addNBR);

var post_fire = l8.merge(l9)
  .filterBounds(aoi)
  .filterDate(post_fireSTART,post_fireEND)
  .map(maskSrClouds)
  .map(addNBR);

// Create varaibles for reference, projection and transformation
var refBand = pre_fire.first().select('SR_B4');
var refProj = refBand.projection();
var refCRS = refProj.crs().getInfo();

var pre_median = pre_fire.median()
var post_median = post_fire.median()

var dNBR = pre_median.select('NBR')
              .subtract(post_median.select('NBR'))
              .rename('dNBR')
              .clip(aoi);

var dem = ee.Image("USGS/SRTMGL1_003");

var hillshade = ee.Terrain.hillshade(dem);

hillshade = hillshade.clip(aoi);

Map.centerObject(aoi, 10);
Map.addLayer(hillshade, {
  min: 0,
  max: 255,
  opacity: 0.35
}, 'Hillshade');


Map.addLayer(dNBR, {
  min: -0.9,
  max: 1.2,
palette: [
'#F7F7F7',
'#C6DBEF',
'#6BAED6',
'#2171B5',
'#08306B'
]
}, 'dNBR');

// Line Fire iamge
var fire_image = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_040036_20240909')
  .clip(aoi);

Map.addLayer(fire_image, {
  bands: ['SR_B7', 'SR_B6', 'SR_B4'],
  min: 7000,
  max: 24000
}, 'Fire (false) Composite');

// Uncomment to export
// Export.image.toDrive({
//   image: dNBR,
//   description: 'dNBR_2024_Fire',
//   folder: exportFolder,
//   fileNamePrefix: 'dNBR_comp',
//   scale: 30,
//   crs: refCRS,
//   region: aoi,
//   maxPixels: 1e13
// });
