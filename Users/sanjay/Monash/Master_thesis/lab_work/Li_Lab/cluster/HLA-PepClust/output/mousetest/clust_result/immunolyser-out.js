
            fetch('corr-data/amino_acids_comparison_with_correlation_1of5_HLA_A0207.json')
                .then(response => {
                    if (!response.ok) {
                        throw new Error("JSON fetch failed");
                    }
                    return response.json();
                })
                .then(data => {
                    console.log("Loaded JSON:", 'corr-data/amino_acids_comparison_with_correlation_1of5_HLA_A0207.json', data);
                    var opt = {"renderer": "canvas", "actions": false};
                    vegaEmbed("#correlation_chart_1of5", data, opt);
                })
                .catch(error => {
                    console.error("Error fetching JSON, loading fallback image instead:", error);
                    document.getElementById("correlation_chart_1of5").innerHTML = `<img src="corr-data/amino_acids_comparison_with_correlation_1of5_HLA_A0207.png" alt="Fallback Image" width="100%">`;
                });
        
            fetch('corr-data/amino_acids_comparison_with_correlation_2of5_HLA_A1102.json')
                .then(response => {
                    if (!response.ok) {
                        throw new Error("JSON fetch failed");
                    }
                    return response.json();
                })
                .then(data => {
                    console.log("Loaded JSON:", 'corr-data/amino_acids_comparison_with_correlation_2of5_HLA_A1102.json', data);
                    var opt = {"renderer": "canvas", "actions": false};
                    vegaEmbed("#correlation_chart_2of5", data, opt);
                })
                .catch(error => {
                    console.error("Error fetching JSON, loading fallback image instead:", error);
                    document.getElementById("correlation_chart_2of5").innerHTML = `<img src="corr-data/amino_acids_comparison_with_correlation_2of5_HLA_A1102.png" alt="Fallback Image" width="100%">`;
                });
        
            fetch('corr-data/amino_acids_comparison_with_correlation_3of5_HLA_C0401.json')
                .then(response => {
                    if (!response.ok) {
                        throw new Error("JSON fetch failed");
                    }
                    return response.json();
                })
                .then(data => {
                    console.log("Loaded JSON:", 'corr-data/amino_acids_comparison_with_correlation_3of5_HLA_C0401.json', data);
                    var opt = {"renderer": "canvas", "actions": false};
                    vegaEmbed("#correlation_chart_3of5", data, opt);
                })
                .catch(error => {
                    console.error("Error fetching JSON, loading fallback image instead:", error);
                    document.getElementById("correlation_chart_3of5").innerHTML = `<img src="corr-data/amino_acids_comparison_with_correlation_3of5_HLA_C0401.png" alt="Fallback Image" width="100%">`;
                });
        
            fetch('corr-data/amino_acids_comparison_with_correlation_4of5_HLA_A0211.json')
                .then(response => {
                    if (!response.ok) {
                        throw new Error("JSON fetch failed");
                    }
                    return response.json();
                })
                .then(data => {
                    console.log("Loaded JSON:", 'corr-data/amino_acids_comparison_with_correlation_4of5_HLA_A0211.json', data);
                    var opt = {"renderer": "canvas", "actions": false};
                    vegaEmbed("#correlation_chart_4of5", data, opt);
                })
                .catch(error => {
                    console.error("Error fetching JSON, loading fallback image instead:", error);
                    document.getElementById("correlation_chart_4of5").innerHTML = `<img src="corr-data/amino_acids_comparison_with_correlation_4of5_HLA_A0211.png" alt="Fallback Image" width="100%">`;
                });
        
            fetch('corr-data/amino_acids_comparison_with_correlation_5of5_HLA_C1701.json')
                .then(response => {
                    if (!response.ok) {
                        throw new Error("JSON fetch failed");
                    }
                    return response.json();
                })
                .then(data => {
                    console.log("Loaded JSON:", 'corr-data/amino_acids_comparison_with_correlation_5of5_HLA_C1701.json', data);
                    var opt = {"renderer": "canvas", "actions": false};
                    vegaEmbed("#correlation_chart_5of5", data, opt);
                })
                .catch(error => {
                    console.error("Error fetching JSON, loading fallback image instead:", error);
                    document.getElementById("correlation_chart_5of5").innerHTML = `<img src="corr-data/amino_acids_comparison_with_correlation_5of5_HLA_C1701.png" alt="Fallback Image" width="100%">`;
                });
        
                // heatmap json or PNG loading
                fetch('corr-data/correlation_heatmap.json')
                .then(response => {
                    if (!response.ok) {
                        throw new Error("JSON fetch failed");
                    }
                    return response.json();
                })
                .then(data => {
                    console.log("Loaded JSON:", 'corr-data/correlation_heatmap.json', data);
                    var opt = {"renderer": "canvas", "actions": false};
                    vegaEmbed("#correlation_heatmap", data, opt);
                })
                .catch(error => {
                    console.error("Error fetching JSON, loading fallback image instead:", error);
                    document.getElementById("correlation_heatmap").innerHTML = `<img src="corr-data/correlation_heatmap.png" alt="Fallback Image" width="100%">`;
                });
        


    // Function to create interactive heatmap
function createHeatmap(containerId, data) {
  // Get container width to make the heatmap responsive
  const containerWidth = document.getElementById(containerId).clientWidth;
  
  // Calculate cell size based on container width and number of columns (HLAs)
  const maxCellSize = 40;
  const minCellSize = 20;
  const calculatedCellSize = Math.max(
    minCellSize, 
    Math.min(maxCellSize, (containerWidth - 150) / data.hlas.length)
  );
  
  // Set up dimensions and margins
  const margin = { top: 50, right: 30, bottom: 120, left: 100 };
  const cellSize = calculatedCellSize;
  const width = Math.min(containerWidth, cellSize * data.hlas.length + margin.left + margin.right);
  const height = cellSize * data.clusters.length + margin.top + margin.bottom;
  
  // Clear any existing SVG
  d3.select(`#${containerId}`).html("");
  
  // Create SVG with responsive width
  const svg = d3.select(`#${containerId}`)
    .append('svg')
    .attr('width', width)
    .attr('height', height)
    .attr('viewBox', `0 0 ${width} ${height}`)
    .attr('preserveAspectRatio', 'xMidYMid meet')
    .append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`);
  
  // Create scales
  const x = d3.scaleBand()
    .domain(data.hlas)
    .range([0, cellSize * data.hlas.length])
    .padding(0.05);
  
  const y = d3.scaleBand()
    .domain(data.clusters)
    .range([0, cellSize * data.clusters.length])
    .padding(0.05);
  
  // Color scale based on selected scheme (YlGnBu by default)
  const colorScale = d3.scaleSequential(d3.interpolateYlGnBu)
    .domain([0, 1]);
  
  // Get current threshold from slider or use default
  const thresholdSlider = document.getElementById('thresholdValue');
  const threshold = thresholdSlider ? parseFloat(thresholdSlider.value) : 0.7;
  
  // Create tooltip
  const tooltip = d3.select('body')
    .append('div')
    .attr('class', 'heatmap-tooltip')
    .style('opacity', 0)
    .style('position', 'absolute')
    .style('padding', '10px')
    .style('background', 'rgba(255, 255, 255, 0.95)')
    .style('border-radius', '4px')
    .style('box-shadow', '0 2px 5px rgba(0, 0, 0, 0.2)')
    .style('pointer-events', 'none')
    .style('z-index', '1000')
    .style('max-width', '220px');
  
  // Add X axis labels
  svg.append('g')
    .selectAll('text')
    .data(data.hlas)
    .enter()
    .append('text')
    .attr('x', d => x(d) + x.bandwidth() / 2)
    .attr('y', -10)
    .attr('text-anchor', 'end')
    .attr('transform', d => `rotate(-45, ${x(d) + x.bandwidth() / 2}, -10)`)
    .text(d => d.replace('HLA_', ''))
    .style('font-size', `${Math.min(cellSize * 0.3, 12)}px`)
    .style('font-weight', '500');
  
  // Add Y axis labels
  svg.append('g')
    .selectAll('text')
    .data(data.clusters)
    .enter()
    .append('text')
    .attr('x', -10)
    .attr('y', d => y(d) + y.bandwidth() / 2)
    .attr('text-anchor', 'end')
    .attr('dominant-baseline', 'middle')
    .text(d => d)
    .style('font-size', `${Math.min(cellSize * 0.3, 12)}px`)
    .style('font-weight', '500');

  // Create the heatmap cells
  svg.selectAll('rect')
    .data(data.heatmapData)
    .enter()
    .append('rect')
    .attr('x', d => x(d.hla))
    .attr('y', d => y(d.cluster))
    .attr('width', x.bandwidth())
    .attr('height', y.bandwidth())
    .style('fill', d => d.correlation === 0 ? '#f0f0f0' : colorScale(d.correlation))
    .style('stroke', '#f5f5f7')
    .style('stroke-width', '1px')
    .style('opacity', d => d.correlation >= threshold ? 1 : 0) // Make cells below threshold disappear
    .on('mouseover', function(event, d) {
      // Only show tooltip and highlight if correlation is above threshold
      if (d.correlation >= threshold) {
        // Highlight the cell
        d3.select(this)
          .style('stroke', '#333')
          .style('stroke-width', '2px');
        
        // Show tooltip
        tooltip.transition()
          .duration(200)
          .style('opacity', 0.9);
        
        tooltip.html(`<strong>Cluster:</strong> ${d.cluster}<br>
                    <strong>HLA:</strong> ${d.hla}<br>
                    <strong>Correlation:</strong> ${d.correlation.toFixed(2)}`)
          .style('left', (event.pageX + 10) + 'px')
          .style('top', (event.pageY - 28) + 'px');
      }
    })
    .on('mouseout', function() {
      // Remove highlight
      d3.select(this)
        .style('stroke', '#f5f5f7')
        .style('stroke-width', '1px');
      
      // Hide tooltip
      tooltip.transition()
        .duration(500)
        .style('opacity', 0);
    });
  
  // Add correlation text to cells with correlation > threshold
  // Adjust font size based on cell size
  const fontSize = Math.min(cellSize * 0.4, 11);
  
  svg.selectAll('text.cell-text')
    .data(data.heatmapData.filter(d => d.correlation >= threshold)) // Only show text for cells above threshold
    .enter()
    .append('text')
    .attr('class', 'cell-text')
    .attr('x', d => x(d.hla) + x.bandwidth() / 2)
    .attr('y', d => y(d.cluster) + y.bandwidth() / 2)
    .attr('text-anchor', 'middle')
    .attr('dominant-baseline', 'middle')
    .text(d => d.correlation.toFixed(2))
    .style('fill', d => d.correlation > 0.5 ? 'white' : 'black')
    .style('font-size', `${fontSize}px`)
    .style('font-weight', 'bold');
}

// Function to create the color scale legend
function createColorLegend(containerId) {
  const containerWidth = document.getElementById(containerId).clientWidth;
  const width = Math.min(300, containerWidth * 0.8);
  const height = 50;
  
  // Clear existing content
  d3.select(`#${containerId}`).html("");
  
  const svg = d3.select(`#${containerId}`)
    .append('svg')
    .attr('width', width)
    .attr('height', height)
    .attr('viewBox', `0 0 ${width} ${height}`)
    .attr('preserveAspectRatio', 'xMidYMid meet');
  
  // Create gradient
  const defs = svg.append('defs');
  const gradient = defs.append('linearGradient')
    .attr('id', 'correlation-gradient')
    .attr('x1', '0%')
    .attr('y1', '0%')
    .attr('x2', '100%')
    .attr('y2', '0%');
  
  // Color stops for YlGnBu
  gradient.append('stop')
    .attr('offset', '0%')
    .attr('stop-color', '#ffffd9');
  
  gradient.append('stop')
    .attr('offset', '20%')
    .attr('stop-color', '#edf8b1');
  
  gradient.append('stop')
    .attr('offset', '40%')
    .attr('stop-color', '#7fcdbb');
  
  gradient.append('stop')
    .attr('offset', '60%')
    .attr('stop-color', '#41b6c4');
  
  gradient.append('stop')
    .attr('offset', '80%')
    .attr('stop-color', '#1d91c0');
  
  gradient.append('stop')
    .attr('offset', '100%')
    .attr('stop-color', '#225ea8');
  
  // Append rectangle filled with gradient
  svg.append('rect')
    .attr('x', 0)
    .attr('y', 10)
    .attr('width', width)
    .attr('height', 20)
    .style('fill', 'url(#correlation-gradient)')
    .style('stroke', '#ccc')
    .style('stroke-width', 1);
  
  // Add labels
  svg.append('text')
    .attr('x', 0)
    .attr('y', 45)
    .text('0.0')
    .style('font-size', '12px')
    .style('text-anchor', 'start');
  
  svg.append('text')
    .attr('x', width / 2)
    .attr('y', 45)
    .text('0.5')
    .style('font-size', '12px')
    .style('text-anchor', 'middle');
  
  svg.append('text')
    .attr('x', width)
    .attr('y', 45)
    .text('1.0')
    .style('font-size', '12px')
    .style('text-anchor', 'end');
}

// Function to handle window resize events
function handleResize(data) {
  if (data) {
    createHeatmap('heatmap-container', data);
    createColorLegend('legend-container');
  }
}

// Function to initialize the heatmap
function initHeatmap() {
  fetch('hla_correlation_data.json')
    .then(response => response.json())
    .then(data => {
      // Store data for resize events
      window.heatmapData = data;
      
      // Create the main heatmap
      createHeatmap('heatmap-container', data);
      
      // Create the color legend
      createColorLegend('legend-container');
      
      // Set up controls
      setupHeatmapControls(data);
      
      // Handle window resize
      window.addEventListener('resize', () => handleResize(data));
    })
    .catch(error => {
      console.error('Error loading the data:', error);
      document.getElementById('heatmap-container').innerHTML = 
        '<div class="alert alert-danger">Error loading heatmap data. See console for details.</div>';
    });
}

// Function to set up interactive controls
function setupHeatmapControls(data) {
  // Color scheme selector
  const colorSchemeSelect = document.getElementById('colorScheme');
  if (colorSchemeSelect) {
    colorSchemeSelect.addEventListener('change', function() {
      const newScheme = this.value;
      updateHeatmapColorScheme(data, newScheme);
    });
  }
  
  // Threshold slider
  const thresholdSlider = document.getElementById('thresholdValue');
  const thresholdDisplay = document.getElementById('thresholdDisplay');
  if (thresholdSlider && thresholdDisplay) {
    thresholdSlider.addEventListener('input', function() {
      const threshold = parseFloat(this.value);
      thresholdDisplay.textContent = threshold.toFixed(2);
      updateHeatmapThreshold(data, threshold);
    });
  }
  
  // Sort selector
  const sortBySelect = document.getElementById('sortBy');
  if (sortBySelect) {
    sortBySelect.addEventListener('change', function() {
      const sortBy = this.value;
      updateHeatmapSorting(data, sortBy);
    });
  }
}

// Function to update color scheme
function updateHeatmapColorScheme(data, scheme) {
  let colorScale;
  
  if (scheme === 'RdYlBu_r') {
    colorScale = d3.scaleSequential(d3.interpolateRdYlBu)
      .domain([1, 0]); // Reversed for this color scheme
  } else if (scheme === 'viridis') {
    colorScale = d3.scaleSequential(d3.interpolateViridis)
      .domain([0, 1]);
  } else if (scheme === 'magma') {
    colorScale = d3.scaleSequential(d3.interpolateMagma)
      .domain([0, 1]);
  } else if (scheme === 'plasma') {
    colorScale = d3.scaleSequential(d3.interpolatePlasma)
      .domain([0, 1]);
  } else {
    // Default YlGnBu
    colorScale = d3.scaleSequential(d3.interpolateYlGnBu)
      .domain([0, 1]);
  }
  
  // Update the cell colors
  d3.selectAll('#heatmap-container rect')
    .style('fill', d => d.correlation === 0 ? '#f0f0f0' : colorScale(d.correlation));
}

// Function to update threshold
function updateHeatmapThreshold(data, threshold) {
  // Make cells below threshold completely disappear (opacity 0)
  d3.selectAll('#heatmap-container rect')
    .style('opacity', d => d.correlation >= threshold ? 1 : 0);
  
  // Redraw the heatmap to update the text labels
  createHeatmap('heatmap-container', data);
}

// Function to update sorting
function updateHeatmapSorting(data, sortBy) {
  // Get current data
  let sortedClusters = [...data.clusters];
  
  if (sortBy === 'correlation') {
    // Calculate average correlation for each cluster
    const clusterAvgCorr = {};
    data.clusters.forEach(cluster => {
      const clusterData = data.heatmapData.filter(item => 
        item.cluster === cluster && item.correlation > 0);
      
      if (clusterData.length > 0) {
        const sum = clusterData.reduce((total, item) => total + item.correlation, 0);
        clusterAvgCorr[cluster] = sum / clusterData.length;
      } else {
        clusterAvgCorr[cluster] = 0;
      }
    });
    
    // Sort by average correlation (descending)
    sortedClusters.sort((a, b) => clusterAvgCorr[b] - clusterAvgCorr[a]);
  } else {
    // Sort by name
    sortedClusters.sort();
  }
  
  // Simply recreate the heatmap with the new sorted clusters
  createHeatmap('heatmap-container', {
    ...data,
    clusters: sortedClusters
  });
}

// Initialize the heatmap when the page loads
document.addEventListener('DOMContentLoaded', initHeatmap);
