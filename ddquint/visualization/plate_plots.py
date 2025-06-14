#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plate plot visualization module for ddQuint with config integration and buffer zone support.

Creates composite plate images showing all 96 wells with individual scatter plots
and appropriate highlighting for aneuploidies and buffer zones. Integrates with
existing clustering results to avoid recomputation.
"""

import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import logging
from matplotlib.ticker import MultipleLocator
from tqdm import tqdm

from ..config import Config, VisualizationError
from ..visualization.well_plots import create_well_plot

logger = logging.getLogger(__name__)


def create_composite_image(results, output_path):
    """
    Create a composite image using existing clustering results without re-running analysis.
    
    Args:
        results (list): List of result dictionaries with clustering data
        output_path (str): Path to save the composite image
        
    Returns:
        str: Path to the saved composite image
        
    Raises:
        VisualizationError: If composite image creation fails
        ValueError: If results data is invalid
    """
    config = Config.get_instance()
    
    if not results:
        error_msg = "No results provided for composite image creation"
        logger.error(error_msg)
        raise ValueError(error_msg)
    
    logger.debug(f"Creating composite image for {len(results)} wells")
    logger.debug(f"Output path: {output_path}")
    
    # Keep track of temporary files
    temp_files = []
    
    try:
        # Get plate layout from config
        row_labels = config.PLATE_ROWS
        col_labels = config.PLATE_COLS
        
        # Generate optimized images for each well
        _generate_well_images(results, output_path, temp_files)
        
        # Create the composite figure
        _create_composite_figure(results, output_path, row_labels, col_labels, config)
        
        logger.debug(f"Composite image saved to: {output_path}")
        return output_path
        
    except Exception as e:
        error_msg = f"Error creating composite image: {str(e)}"
        logger.error(error_msg)
        logger.debug(f"Error details: {str(e)}", exc_info=True)
        raise VisualizationError(error_msg) from e
        
    finally:
        # Clean up temporary files
        _cleanup_temp_files(temp_files, results)


def _generate_well_images(results, output_path, temp_files):
    """Generate optimized images for each well with progress bar."""
    config = Config.get_instance()
    
    for result in tqdm(results, desc="Creating Plate image", unit="well"):
        if not result.get('well'):
            continue
            
        # Get the data file path
        df_file = _get_data_file_path(result, config)
        if not df_file or not os.path.exists(df_file):
            logger.debug(f"Raw data file not found: {df_file}")
            continue
            
        # Load and process the raw data
        df_clean = _load_and_clean_data(df_file)
        if df_clean is None:
            continue
            
        # Use existing clustering results
        clustering_results = _extract_clustering_results(result)
        if not _validate_clustering_results(clustering_results, result['well']):
            continue
        
        # Create optimized plot for composite image
        temp_path = _create_temp_well_plot(result, output_path, df_clean, 
                                         clustering_results, config, temp_files)
        if temp_path:
            result['temp_graph_path'] = temp_path


def _get_data_file_path(result, config):
    """Get the path to the raw data file for a result."""
    try:
        return os.path.join(os.path.dirname(result['graph_path']), "..", 
                           config.RAW_DATA_DIR_NAME, result['filename'])
    except (KeyError, TypeError):
        logger.debug(f"Could not construct data file path for result: {result.get('well', 'unknown')}")
        return None


def _load_and_clean_data(df_file):
    """Load and clean CSV data file."""
    try:
        # Find the header row
        header_row = None
        with open(df_file, 'r', encoding='utf-8', errors='ignore') as f:
            for i, line in enumerate(f):
                if ('Ch1Amplitude' in line or 'Ch1 Amplitude' in line) and \
                   ('Ch2Amplitude' in line or 'Ch2 Amplitude' in line):
                    header_row = i
                    break
        
        if header_row is None:
            logger.debug(f"Could not find header row in {os.path.basename(df_file)}")
            return None
            
        # Load the CSV data
        df = pd.read_csv(df_file, skiprows=header_row)
        
        # Check for required columns
        required_cols = ['Ch1Amplitude', 'Ch2Amplitude']
        if not all(col in df.columns for col in required_cols):
            logger.debug(f"Required columns not found in {os.path.basename(df_file)}")
            return None
        
        # Filter rows with NaN values and create explicit copy
        df_clean = df[required_cols].dropna().copy()
        return df_clean
        
    except Exception as e:
        logger.debug(f"Error loading data file {os.path.basename(df_file)}: {e}")
        return None


def _extract_clustering_results(result):
    """Extract clustering results from result dictionary."""
    return {
        'df_filtered': result.get('df_filtered'),
        'target_mapping': result.get('target_mapping'),
        'counts': result.get('counts', {}),
        'copy_numbers': result.get('copy_numbers', {}),
        'copy_number_states': result.get('copy_number_states', {}),
        'has_aneuploidy': result.get('has_aneuploidy', False),
        'has_buffer_zone': result.get('has_buffer_zone', False),
        'chrom3_reclustered': result.get('chrom3_reclustered', False)
    }


def _validate_clustering_results(clustering_results, well_id):
    """Validate that clustering results contain necessary data."""
    if clustering_results['df_filtered'] is None or clustering_results['target_mapping'] is None:
        logger.debug(f"Missing clustering data for well {well_id}")
        return False
    return True


def _create_temp_well_plot(result, output_path, df_clean, clustering_results, config, temp_files):
    """Create temporary well plot for composite image."""
    try:
        output_dir = os.path.dirname(output_path)
        graphs_dir = os.path.join(output_dir, config.GRAPHS_DIR_NAME)
        os.makedirs(graphs_dir, exist_ok=True)
        
        # Create temp file in the Graphs directory
        temp_path = os.path.join(graphs_dir, f"{result['well']}_temp.png")
        create_well_plot(df_clean, clustering_results, result['well'], 
                         temp_path, for_composite=True, add_copy_numbers=True)
        
        # Track the temporary file
        temp_files.append(temp_path)
        return temp_path
        
    except Exception as e:
        logger.debug(f"Error creating temp plot for well {result.get('well')}: {e}")
        return None


def _create_composite_figure(results, output_path, row_labels, col_labels, config):
    """Create the main composite figure."""
    # Create figure with configured size
    fig_size = config.COMPOSITE_FIGURE_SIZE
    fig = plt.figure(figsize=fig_size)
    logger.debug(f"Creating composite figure with size: {fig_size}")
    
    # Create GridSpec with spacing
    gs = gridspec.GridSpec(8, 12, figure=fig, wspace=0.02, hspace=0.02)
    
    # Create well results mapping
    well_results = {r['well']: r for r in results if r.get('well') is not None}
    
    # Ensure proper margins
    plt.subplots_adjust(left=0.04, right=0.96, top=0.96, bottom=0.04)
    
    # Create subplot for each well position
    _create_well_subplots(fig, gs, row_labels, col_labels, well_results, config)
    
    # Add row and column labels
    _add_plate_labels(fig, row_labels, col_labels)
    
    # Save the composite image
    fig.savefig(output_path, dpi=400, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)


def _create_well_subplots(fig, gs, row_labels, col_labels, well_results, config):
    """Create individual subplots for each well position."""
    for i, row in enumerate(row_labels):
        for j, col_num in enumerate(range(1, int(col_labels[-1]) + 1)):
            col = str(col_num)
            well = config.WELL_FORMAT.format(row=row, col=int(col))
            
            # Add subplot at this position
            ax = fig.add_subplot(gs[i, j])
            ax.set_facecolor('#f5f5f5')  # Light gray background
            
            if well in well_results:
                _populate_data_well(ax, well, well_results[well], config)
            else:
                _populate_empty_well(ax, well, config)
            
            # Keep axis visibility for all plots
            ax.set_xticks([])
            ax.set_yticks([])


def _populate_data_well(ax, well, result, config):
    """Populate subplot for a well with data."""
    # Use temp graph path if available, otherwise fall back to original
    graph_path = result.get('temp_graph_path', result.get('graph_path'))
    
    if graph_path and os.path.exists(graph_path):
        try:
            # Read and display the individual well image
            img = plt.imread(graph_path)
            ax.imshow(img)
            
            # Add title
            sample_name = result.get('sample_name')
            title = sample_name if sample_name else well
            ax.set_title(title, fontsize=6, pad=2)
            
            # Apply colored borders based on copy number state
            _apply_well_border(ax, result, well)
            
        except Exception as e:
            logger.debug(f"Error displaying image for well {well}: {e}")
            ax.text(0.5, 0.5, "Image Error", 
                    horizontalalignment='center', verticalalignment='center', 
                    transform=ax.transAxes, fontsize=8, color='red')
    else:
        ax.text(0.5, 0.5, "No Image", 
                horizontalalignment='center', verticalalignment='center', 
                transform=ax.transAxes, fontsize=8)


def _populate_empty_well(ax, well, config):
    """Populate subplot for an empty well."""
    # Create properly sized placeholder
    axis_limits = config.get_axis_limits()
    ax.set_xlim(axis_limits['x'])
    ax.set_ylim(axis_limits['y'])
    
    # Add grid with configured spacing
    ax.grid(True, alpha=0.4, linewidth=0.8)
    grid_intervals = config.get_grid_intervals()
    ax.xaxis.set_major_locator(MultipleLocator(grid_intervals['x']))
    ax.yaxis.set_major_locator(MultipleLocator(grid_intervals['y']))
    
    # Turn off tick marks but keep the grid
    ax.tick_params(axis='both', which='both', length=0)
    ax.set_aspect('auto')
    
    # Add well identifier in gray
    ax.text(0.5, 0.5, well, fontsize=8, color='gray',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    
    # Apply default grey border
    for spine in ax.spines.values():
        spine.set_color('#cccccc')
        spine.set_linewidth(1)
        spine.set_visible(True)


def _apply_well_border(ax, result, well):
    """Apply colored border based on well status."""
    # Buffer zone trumps aneuploidy in detection
    if result.get('has_buffer_zone', False):
        border_color = '#000000'  # Black border for buffer zone
        border_width = 1
        logger.debug(f"Applied buffer zone border (black) to well {well}")
    elif result.get('has_aneuploidy', False):
        border_color = '#E6B8E6'  # Light purple border for aneuploidy
        border_width = 3
        logger.debug(f"Applied aneuploidy border (light purple) to well {well}")
    else:
        border_color = '#B0B0B0'  # Light grey border for euploid
        border_width = 1
        logger.debug(f"Applied euploid border (light grey) to well {well}")
    
    # Apply the border
    for spine in ax.spines.values():
        spine.set_edgecolor(border_color)
        spine.set_color(border_color)
        spine.set_linewidth(border_width)
        spine.set_visible(True)


def _add_plate_labels(fig, row_labels, col_labels):
    """Add row and column labels to the plate."""
    # Add row labels (A-H) with proper alignment
    for i, row in enumerate(row_labels):
        ax = fig.axes[i * 12]  # Get the first plot in this row
        y_center = (ax.get_position().y0 + ax.get_position().y1) / 2
        fig.text(0.02, y_center, row, ha='center', va='center', fontsize=12, weight='bold')
    
    # Add column labels (1-12) with proper alignment
    for j, col in enumerate(col_labels):
        ax = fig.axes[j]  # Get the plot in the first row for this column
        x_center = (ax.get_position().x0 + ax.get_position().x1) / 2
        fig.text(x_center, 0.98, col, ha='center', va='center', fontsize=12, weight='bold')


def _cleanup_temp_files(temp_files, results):
    """Clean up temporary files and references."""
    for temp_file in temp_files:
        try:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                logger.debug(f"Removed temporary file: {os.path.basename(temp_file)}")
        except Exception as e:
            logger.debug(f"Error deleting temporary file {os.path.basename(temp_file)}: {e}")
    
    # Clear any references to temporary files in results
    for result in results:
        if 'temp_graph_path' in result:
            del result['temp_graph_path']