#!/usr/bin/env python3
from typing import List, Dict, Optional, Tuple, Set
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import argparse
from pathlib import Path
import logging
import gzip
from collections import defaultdict
import numpy as np
import re
import json
import plotly.express as px

def parse_gff_attributes(attr_str: str) -> Dict[str, str]:
    """Parse GFF3 attributes field"""
    if attr_str == '.':
        return {}
    
    attributes = {}
    for attr in attr_str.strip(';').split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key.strip()] = value.strip()
    return attributes

def parse_gff_line(line: str) -> Optional[Dict]:
    """Parse a single GFF3 line"""
    if line.startswith('#') or not line.strip():
        return None
        
    fields = line.strip().split('\t')
    if len(fields) != 9:
        return None
        
    try:
        return {
            'seqid': fields[0],
            'source': fields[1],
            'type': fields[2],
            'start': int(fields[3]),
            'end': int(fields[4]),
            'score': fields[5],
            'strand': fields[6],
            'phase': fields[7],
            'attributes': parse_gff_attributes(fields[8])
        }
    except (ValueError, IndexError):
        return None

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class GenomeVisualizer:
    """Class for visualizing genome structure from GFF files"""
    
    FEATURE_COLORS = {
        'telomere': '#FF9999',
        'subtelomere': '#FFB366',
        'centromere': '#FF99CC',
        'pericentromere': '#FF99FF',
        'chromatin_region': '#99FF99',
        'TAD': '#99FFFF',
        'satellite': '#9999FF',
        'mobile_element': '#FFFF99',
        'CTCF_site': '#FF99FF',
        'chromatin_transition': '#FFCC99'
    }
    
    def __init__(self, gff_path: Path):
        """Initialize visualizer with GFF file path"""
        self.gff_path = gff_path
        self.features = self._load_gff()
        self.chromosome_length = max(self.features['end'])
        
    def _load_gff(self) -> pd.DataFrame:
        """Load and parse GFF file into DataFrame"""
        # Open GFF file (handle both compressed and uncompressed)
        if str(self.gff_path).endswith('.gz'):
            open_func = gzip.open
        else:
            open_func = open
        
        features = []
        with open_func(self.gff_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                rec = parse_gff_line(line)
                if rec:
                    feat = {
                        'type': rec['type'],
                        'start': rec['start'],
                        'end': rec['end'],
                        'strand': rec['strand'],
                        'score': rec['score']
                    }
                    # Add all attributes
                    feat.update(rec['attributes'])
                    
                    # Convert numeric attributes
                    for key in ['gc_content', 'gradient_value', 'strength']:
                        if key in feat:
                            try:
                                feat[key] = float(feat[key])
                            except ValueError:
                                pass
                                
                    features.append(feat)
        
        return pd.DataFrame(features)
    
    def _create_ideogram(self) -> go.Figure:
        """Create chromosome ideogram"""
        fig = go.Figure()
        
        # Add features by type
        for feature_type in self.features['type'].unique():
            type_features = self.features[self.features['type'] == feature_type]
            
            # Get color for feature type
            color = self.FEATURE_COLORS.get(feature_type, '#CCCCCC')
            
            # Add rectangles for each feature
            fig.add_trace(go.Scatter(
                x=type_features['start'].tolist() + type_features['end'].tolist()[::-1],
                y=[feature_type] * len(type_features) * 2,
                fill='toself',
                fillcolor=color,
                line=dict(color='black', width=1),
                name=feature_type,
                showlegend=True
            ))
        
        # Update layout
        fig.update_layout(
            title='Chromosome Structure Ideogram',
            xaxis_title='Position (bp)',
            yaxis_title='Feature Type',
            height=400,
            showlegend=True,
            plot_bgcolor='white',
            yaxis=dict(
                showgrid=False,
                zeroline=False
            ),
            xaxis=dict(
                showgrid=True,
                gridwidth=1,
                gridcolor='lightgray',
                zeroline=False
            )
        )
        
        return fig
    
    def _create_feature_density(self) -> go.Figure:
        """Create feature density plot"""
        # Create bins for density calculation
        bin_size = self.chromosome_length // 1000  # 1000 bins
        bins = np.linspace(0, self.chromosome_length, 1001)
        
        # Create subplot figure
        fig = make_subplots(rows=len(self.features['type'].unique()), cols=1,
                           shared_xaxes=True,
                           subplot_titles=self.features['type'].unique())
        
        # Plot density for each feature type
        for i, feature_type in enumerate(self.features['type'].unique(), 1):
            type_features = self.features[self.features['type'] == feature_type]
            
            # Calculate feature density
            densities, _ = np.histogram(
                np.concatenate([
                    np.linspace(start, end, 100) 
                    for start, end in zip(type_features['start'], type_features['end'])
                ]),
                bins=bins
            )
            
            # Add density plot
            fig.add_trace(
                go.Scatter(
                    x=(bins[:-1] + bins[1:]) / 2,
                    y=densities,
                    fill='tozeroy',
                    name=feature_type,
                    fillcolor=self.FEATURE_COLORS.get(feature_type, '#CCCCCC'),
                    line=dict(color='black', width=1)
                ),
                row=i, col=1
            )
        
        # Update layout
        fig.update_layout(
            title='Feature Density Distribution',
            height=200 * len(self.features['type'].unique()),
            showlegend=False,
            plot_bgcolor='white'
        )
        
        # Update all xaxes
        fig.update_xaxes(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=False,
            title_text='Position (bp)'
        )
        
        # Update all yaxes
        fig.update_yaxes(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=False,
            title_text='Density'
        )
        
        return fig
    
    def _create_gc_content_plot(self) -> go.Figure:
        """Create GC content distribution plot"""
        gc_features = self.features[self.features['gc_content'].notna()]
        
        if len(gc_features) == 0:
            return None
        
        fig = go.Figure()
        
        # Add scatter plot for GC content
        for feature_type in gc_features['type'].unique():
            type_features = gc_features[gc_features['type'] == feature_type]
            
            fig.add_trace(go.Scatter(
                x=(type_features['start'] + type_features['end']) / 2,
                y=type_features['gc_content'],
                mode='markers',
                name=feature_type,
                marker=dict(
                    color=self.FEATURE_COLORS.get(feature_type, '#CCCCCC'),
                    size=8
                )
            ))
        
        # Update layout
        fig.update_layout(
            title='GC Content Distribution',
            xaxis_title='Position (bp)',
            yaxis_title='GC Content',
            height=400,
            showlegend=True,
            plot_bgcolor='white',
            yaxis=dict(
                showgrid=True,
                gridwidth=1,
                gridcolor='lightgray',
                zeroline=False,
                range=[0, 1]
            ),
            xaxis=dict(
                showgrid=True,
                gridwidth=1,
                gridcolor='lightgray',
                zeroline=False
            )
        )
        
        return fig
    
    def _create_feature_statistics(self) -> go.Figure:
        """Create feature statistics visualization"""
        # Calculate basic statistics
        stats = []
        for feature_type in self.features['type'].unique():
            type_features = self.features[self.features['type'] == feature_type]
            
            stats.append({
                'type': feature_type,
                'count': len(type_features),
                'total_length': sum(type_features['end'] - type_features['start']),
                'avg_length': np.mean(type_features['end'] - type_features['start']),
                'coverage': sum(type_features['end'] - type_features['start']) / self.chromosome_length
            })
        
        stats_df = pd.DataFrame(stats)
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Feature Counts', 'Total Length', 'Average Length', 'Chromosome Coverage')
        )
        
        # Add bar plots
        # Feature counts
        fig.add_trace(
            go.Bar(
                x=stats_df['type'],
                y=stats_df['count'],
                marker_color=[self.FEATURE_COLORS.get(t, '#CCCCCC') for t in stats_df['type']]
            ),
            row=1, col=1
        )
        
        # Total length
        fig.add_trace(
            go.Bar(
                x=stats_df['type'],
                y=stats_df['total_length'],
                marker_color=[self.FEATURE_COLORS.get(t, '#CCCCCC') for t in stats_df['type']]
            ),
            row=1, col=2
        )
        
        # Average length
        fig.add_trace(
            go.Bar(
                x=stats_df['type'],
                y=stats_df['avg_length'],
                marker_color=[self.FEATURE_COLORS.get(t, '#CCCCCC') for t in stats_df['type']]
            ),
            row=2, col=1
        )
        
        # Coverage
        fig.add_trace(
            go.Bar(
                x=stats_df['type'],
                y=stats_df['coverage'],
                marker_color=[self.FEATURE_COLORS.get(t, '#CCCCCC') for t in stats_df['type']]
            ),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            height=800,
            title='Feature Statistics',
            showlegend=False,
            plot_bgcolor='white'
        )
        
        # Update all xaxes
        fig.update_xaxes(
            tickangle=45,
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=False
        )
        
        # Update all yaxes
        fig.update_yaxes(
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=False
        )
        
        return fig
    
    def generate_report(self, output_dir: Path):
        """Generate comprehensive visualization report"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create visualizations
        ideogram = self._create_ideogram()
        density = self._create_feature_density()
        gc_content = self._create_gc_content_plot()
        statistics = self._create_feature_statistics()
        
        # Save figures
        ideogram.write_html(output_dir / 'ideogram.html')
        density.write_html(output_dir / 'density.html')
        if gc_content is not None:
            gc_content.write_html(output_dir / 'gc_content.html')
        statistics.write_html(output_dir / 'statistics.html')
        
        # Create index.html
        index_html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Genome Structure Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .plot {{ margin: 20px 0; }}
                iframe {{ border: none; width: 100%; }}
            </style>
        </head>
        <body>
            <h1>Genome Structure Report</h1>
            
            <div class="plot">
                <h2>Chromosome Ideogram</h2>
                <iframe src="ideogram.html" height="450px"></iframe>
            </div>
            
            <div class="plot">
                <h2>Feature Density Distribution</h2>
                <iframe src="density.html" height="{250 * len(self.features['type'].unique())}px"></iframe>
            </div>
            
            {'<div class="plot"><h2>GC Content Distribution</h2><iframe src="gc_content.html" height="450px"></iframe></div>' 
             if gc_content is not None else ''}
            
            <div class="plot">
                <h2>Feature Statistics</h2>
                <iframe src="statistics.html" height="850px"></iframe>
            </div>
        </body>
        </html>
        """
        
        with open(output_dir / 'index.html', 'w') as f:
            f.write(index_html)
        
        logger.info(f"Generated visualization report in {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Visualize genome structure from GFF file')
    parser.add_argument('gff_file', type=Path, help='Path to GFF file')
    parser.add_argument('--output-dir', type=Path, default=Path('viz_output'),
                        help='Output directory (default: viz_output)')
    
    args = parser.parse_args()
    
    try:
        visualizer = GenomeVisualizer(args.gff_file)
        visualizer.generate_report(args.output_dir)
        
    except Exception as e:
        logger.error(f"Error during visualization: {e}", exc_info=True)
        raise

if __name__ == '__main__':
    main()