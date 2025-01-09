# LST_Qena
calculating the LST for the archaeological Sheikh Al Arab Hamam area
//vis params
var vizParams = {
  bands: ['B5', 'B6', 'B4'],
  min: 642,
  max: 3307,
  gamma: [1, 0.9, 1.1]
};

var vizParams2 = {
  bands: ['B4', 'B3', 'B2'],
  min: 0,
  max: 3000,
  gamma: 1.4,
};

//RADIANCE_MULT_BAND_6_VCID_2: 0.03720499947667122
//RADIANCE_ADD_BAND_6_VCID_2: 3.1628000736236572
//load the collection:
{
var col = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA').filterMetadata ('CLOUD_COVER', 'less_than', 5)
    .filterDate('2013-01-01','2013-12-31')
    .filterBounds(geometry)
    .map(function(image){return image.clip(geometry)});
}


print('collection', col);

Map.addLayer(col.first());

Map.centerObject(geometry);
//imagen reduction

var image = col.median();
//print('image', image);

//Map.addLayer(image, vizParams2);

//median
var ndvi1 = image.normalizedDifference(['B5', 'B4']).rename('NDVI');
var ndviParams = {min: 0.10554729676864096, max: 0.41295681063122924, palette: ['blue', 'white', 'green']};

//print('ndvi1', ndvi1);

//individual LST images

var col_list = col.toList(col.size());

var LST_col = col_list.map(function (ele) {
  
  var date = ee.Image(ele).get('system:time_start');

  var ndvi = ee.Image(ele).normalizedDifference(['B5', 'B4']).rename('NDVI');
  
  // find the min and max of NDVI
  var min = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.min(),
    geometry: geometry,
    scale: 10,
    maxPixels: 1e9,
    bestEffort:true
  }).values().get(0));
  
  var max = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: geometry,
    scale: 10,
    maxPixels: 1e9
  }).values().get(0));
  
  var fv = (ndvi.subtract(min).divide(max.subtract(min))).pow(ee.Number(2)).rename('FV');
  
  var a= ee.Number(0.004);
  var b= ee.Number(0.986);
  
  var EM = fv.multiply(a).add(b).rename('EMM');

  var image = ee.Image(ele);

  var Tb = image.expression(
    '1282.7099609375 / log ((666.0900268554688/(0.03720499947667122*L + 3.1628000736236572)) + 1 )', {
    'L': image.select('B10')
  });

  var LST = image.expression(
    '(Tb/(1 + (0.00115* (Tb / 1.438))*log(Ep)))-273.15', {
      'Tb': Tb,
      'Ep': fv.multiply(a).add(b)
  });

  return ee.Algorithms.If(min, LST.set('system:time_start', date).float().rename('LST'), 0);

}).removeAll([0]);

LST_col = ee.ImageCollection(LST_col);

print("LST_col", LST_col);

// Calculate Normalized Difference Vegetation Index (NDVI)
var ndvi = image.normalizedDifference(['B5', 'B4']).rename('NDVI');

// Define NDVI Visualization Parameters
var ndviPalette = {
 min: -1,
 max: 1,
 palette: ['blue', 'white', 'green']
};

Map.addLayer(ndvi, ndviPalette, 'NDVI Yogyakarta');
// Calculate the minimum NDVI value within the AOI
var ndviMin = ee.Number(ndvi.reduceRegion({
  reducer   : ee.Reducer.min(),
  geometry  : geometry,
  scale     : 30,
  maxPixels : 1e9
}).values().get(0));

// Calculate the maximum NDVI value within the AOI
var ndviMax = ee.Number(ndvi.reduceRegion({
  reducer   : ee.Reducer.max(),
  geometry  : geometry,
  scale     : 30,
  maxPixels : 1e9
}).values().get(0));
// Print the Minimum and Maximum NDVI Values
print("Minimum NDVI:", ndviMin);
print("Maximum NDVI:", ndviMax);
// Fraction of Vegetation (FV) Calculation
// Formula: ((NDVI - NDVI_min) / (NDVI_max - NDVI_min))^2

// Calculate the Fraction of Vegetation (FV) using the NDVI values within the specified range.
// NDVI_min represents the minimum NDVI value, and NDVI_max represents the maximum NDVI value
var fv = ((ndvi.subtract(ndviMin)).divide(ndviMax.subtract(ndviMin)))
          .pow(ee.Number(2))
          .rename('FV');
        
  
// Emissivity Calculation
// Formula: 0.004 * FV + 0.986

// Calculate Land Surface Emissivity (EM) using the Fraction of Vegetation (FV).
// The 0.004 coefficient represents the emissivity variation due to vegetation,
// and the 0.986 represents the base emissivity for other surfaces.

var em = fv.multiply(ee.Number(0.004)).add(ee.Number(0.986)).rename('EM');
// Select Thermal Band (Band 10) and Rename It
var thermal = image.select('B10').rename('thermal');
// Now, lets calculate the land surface temperature (LST)

// Formula: (TB / (1 + (λ * (TB / 1.438)) * ln(em))) - 273.15
var lst = thermal.expression(
  '(TB / (1 + (0.00115 * (TB / 1.438)) * log(em))) - 273.15', {
    'TB': thermal.select('thermal'), // Select the thermal band (TB)
    'em': em // Assign emissivity (em)
  }).rename('LST Yogyakarta 2022');

// Add the LST Layer to the Map with Custom Visualization Parameters
Map.addLayer(lst, {
  min: 18.47, // Minimum LST value
  max: 42.86, // Maximum LST value
  palette: [
    '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
    '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
    '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
    'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
    'ff0000', 'de0101', 'c21301', 'a71001', '911003'
  ]}, 'Land Surface Temperature 2022');
  // Create a Legend for Land Surface Temperature (LST) Visualization
var minLST = 15; // Minimum LST value
var maxLST = 45; // Maximum LST value

// Create a panel for the legend with styling
var legend = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px 15px',
    backgroundColor: 'white'
  }
});

// Create a title for the legend
var legendTitle = ui.Label({
  value: 'Land Surface Temperature (°C)',
  style: {
    fontWeight: 'bold',
    fontSize: '20px',
    margin: '0 0 4px 0',
    padding: '0'
  }
});

// Add the legend title to the legend panel
legend.add(legendTitle);

// Define a color palette for the legend
var palette = [
  '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
  '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
  '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
  'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
  'ff0000', 'de0101', 'c21301', 'a71001', '911003', '210300'
];

// Calculate the step value for the legend
var step = (maxLST - minLST) / (palette.length - 1);

// Loop through the palette and create legend entries
for (var i = 0; i < palette.length; i++) {
  // Create a color box for each legend entry
  var colorBox = ui.Label({
    style: {
      backgroundColor: '#' + palette[i],
      padding: '8px',
      margin: '0 0 8px 0',
      width: '42px'
    }
  });

  // Create a label with LST values for each legend entry
  var legendLabel = ui.Label({
    value: (minLST + i * step).toFixed(2),
    style: { margin: '0 0 10px 6px' }
  });

  // Create a panel to arrange color box and label horizontally
  var legendPanel = ui.Panel({
    widgets: [colorBox, legendLabel],
    layout: ui.Panel.Layout.Flow('horizontal')
  });

  // Add the legend entry panel to the main legend panel
  legend.add(legendPanel);
}

// Add the legend to the Google Earth Engine map
Map.add(legend);

// Create a Map Title
var mapTitle = ui.Panel({
  style: {
    position: 'top-center',
    padding: '20px 20px'
  }
});
var mapTitle2 = ui.Label({
  value: 'Land Surface Temperature - 2022 - Yogyakarta',
  style: {
    fontWeight: 'bold',
    fontSize: '30px',
    margin: '8px 8px 8px 8px',
    padding: '0'
  }
});
mapTitle.add(mapTitle2);

// Export the LST image to your Google Drive
Export.image.toDrive({
  image: lst,
  description: 'lst_01-12-2013',
  folder: 'EarthEngineExports',
  scale: 10,
  region: geometry,
  maxPixels: 1e9,
  crs: 'EPSG:4326'
});


