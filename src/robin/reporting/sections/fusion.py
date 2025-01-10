"""
fusion.py

This module contains the fusion analysis section of the report.
"""

import os
import logging
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.units import inch
from .base import ReportSection
import pandas as pd
import networkx as nx
from reportlab.lib.styles import ParagraphStyle

logger = logging.getLogger(__name__)

class FusionSection(ReportSection):
    """Section containing the fusion analysis results."""

    def _get_gene_pairs(self, data):
        """Get unique gene fusion pairs from the data."""
        if data is None or data.empty:
            return []
        
        # Group reads by readID to find pairs
        read_groups = data.groupby('readID')
        gene_pairs = []
        
        # First collect all read IDs for each gene pair
        gene_pair_reads = {}
        for _, group in read_groups:
            genes = group['Gene'].unique()
            if len(genes) >= 2:
                # Sort genes to ensure consistent ordering
                genes = sorted(genes)
                for i in range(len(genes)-1):
                    for j in range(i+1, len(genes)):
                        pair = (genes[i], genes[j])
                        if pair not in gene_pair_reads:
                            gene_pair_reads[pair] = set()
                        gene_pair_reads[pair].add(group['readID'].iloc[0])
        
        # Only include pairs with 3 or more supporting reads
        for pair, reads in gene_pair_reads.items():
            if len(reads) >= 3:
                gene_pairs.append(pair)
        
        # Get unique pairs
        unique_pairs = list(set(gene_pairs))
        return unique_pairs

    def _get_gene_networks(self, gene_pairs):
        """Create networks of connected genes."""
        G = nx.Graph()
        for pair in gene_pairs:
            G.add_edge(pair[0], pair[1])
        connected_components = list(nx.connected_components(G))
        return [list(component) for component in connected_components]

    def _load_fusion_data(self):
        """Load fusion data from the output directory."""
        fusion_data = {
            'master_candidates': None,
            'all_candidates': None,
            'master_path': os.path.join(self.report.output, 'fusion_candidates_master.csv'),
            'all_path': os.path.join(self.report.output, 'fusion_candidates_all.csv')
        }
        
        try:
            if os.path.exists(fusion_data['master_path']):
                fusion_data['master_candidates'] = pd.read_csv(
                    fusion_data['master_path'],
                    names=['chromBED', 'BS', 'BE', 'Gene', 'chrom', 'mS', 'mE', 'readID', 'mapQ', 'strand',
                          'Read_Map_Start', 'Read_Map_End', 'Secondary', 'Supplementary', 'mapping_span']
                )
                logger.debug(f"Loaded master fusion candidates from {fusion_data['master_path']}")
            
            if os.path.exists(fusion_data['all_path']):
                fusion_data['all_candidates'] = pd.read_csv(
                    fusion_data['all_path'],
                    names=['chromBED', 'BS', 'BE', 'Gene', 'chrom', 'mS', 'mE', 'readID', 'mapQ', 'strand',
                          'Read_Map_Start', 'Read_Map_End', 'Secondary', 'Supplementary', 'mapping_span']
                )
                logger.debug(f"Loaded all fusion candidates from {fusion_data['all_path']}")
        
        except Exception as e:
            logger.error(f"Error loading fusion data: {str(e)}")
            logger.debug("Exception details:", exc_info=True)
        
        return fusion_data

    def _format_fusion_table(self, data, title):
        """Create a formatted table for fusion data."""
        if data is None or data.empty:
            return None

        # Group by readID to get fusion pairs
        read_groups = data.groupby('readID')
        
        # Create paragraph styles for cells
        header_style = ParagraphStyle(
            'HeaderStyle',
            parent=self.styles.styles["Normal"],
            fontName='Helvetica-Bold',
            fontSize=8,
            textColor=self.styles.COLORS["primary"],
            leading=10,
            spaceBefore=0,
            spaceAfter=0
        )
        
        cell_style = ParagraphStyle(
            'CellStyle',
            parent=self.styles.styles["Normal"],
            fontName='Helvetica',
            fontSize=8,
            textColor=colors.black,
            leading=10,
            spaceBefore=0,
            spaceAfter=0
        )
        
        # Prepare table data with updated headers using Paragraph objects
        table_data = [[
            Paragraph('Fusion Pair', header_style),
            Paragraph('Chr 1', header_style),
            Paragraph('Chr 2', header_style),
            Paragraph('Gene 1 Position', header_style),
            Paragraph('Gene 2 Position', header_style),
            Paragraph('Reads', header_style)
        ]]
        
        # Track processed pairs to avoid duplicates
        processed_pairs = set()
        fusion_details = {}
        
        # First pass: collect all read IDs for each gene pair
        gene_pair_reads = {}
        for read_id, group in read_groups:
            if len(group) >= 2:
                genes = sorted(group['Gene'].unique())  # Sort to ensure consistent ordering
                if len(genes) >= 2:
                    for i in range(len(genes)-1):
                        for j in range(i+1, len(genes)):
                            gene_pair = f"{genes[i]}-{genes[j]}"
                            if gene_pair not in gene_pair_reads:
                                gene_pair_reads[gene_pair] = set()
                            gene_pair_reads[gene_pair].add(read_id)
        
        # Second pass: process each gene pair with 3 or more supporting reads
        for gene_pair, reads in gene_pair_reads.items():
            if len(reads) >= 3 and gene_pair not in processed_pairs:
                processed_pairs.add(gene_pair)
                genes = gene_pair.split('-')
                
                # Get details for each gene in the pair
                gene_info = {}
                for gene in genes:
                    gene_data = data[data['Gene'] == gene].iloc[0]
                    gene_info[gene] = {
                        'chrom': gene_data['chromBED'],
                        'position': f"{gene_data['BS']}-{gene_data['BE']}"
                    }
                
                # Store fusion details
                fusion_details[gene_pair] = {
                    'chrom1': gene_info[genes[0]]['chrom'],
                    'chrom2': gene_info[genes[1]]['chrom'],
                    'pos1': gene_info[genes[0]]['position'],
                    'pos2': gene_info[genes[1]]['position'],
                    'supporting_reads': len(reads)
                }
        
        # Add rows to table using Paragraph objects for wrappable text
        for gene_pair, details in fusion_details.items():
            table_data.append([
                Paragraph(gene_pair, cell_style),
                Paragraph(details['chrom1'], cell_style),
                Paragraph(details['chrom2'], cell_style),
                Paragraph(details['pos1'], cell_style),
                Paragraph(details['pos2'], cell_style),
                Paragraph(str(details['supporting_reads']), cell_style)
            ])

        # Return None if no fusions meet the threshold
        if len(table_data) == 1:  # Only header row
            return None

        # Define column widths (in inches) - adjusted for better proportions
        col_widths = [1.8, 0.5, 0.5, 1.6, 1.6, 0.6]
        
        # Create and style the table
        table = Table(table_data, colWidths=[inch * width for width in col_widths], repeatRows=1)
        table.setStyle(TableStyle([
            # Inherit modern table style
            *self.MODERN_TABLE_STYLE._cmds,
            
            # Preserve specific alignments
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),    # Left-align Fusion Pair
            ('ALIGN', (1, 0), (2, -1), 'CENTER'),  # Center-align Chromosomes
            ('ALIGN', (3, 0), (4, -1), 'LEFT'),    # Left-align positions
            ('ALIGN', (5, 0), (5, -1), 'RIGHT'),   # Right-align Supporting Reads
        ]))
        
        return table

    def add_content(self):
        """Add the fusion analysis content to the report."""
        logger.debug("Starting fusion section content generation")
        
        # Add section title
        self.elements.append(Paragraph("Gene Fusion Analysis", self.styles.styles["Heading1"]))
        self.elements.append(Spacer(1, 12))
        
        # Load fusion data
        fusion_data = self._load_fusion_data()
        
        # Get unique fusion pairs
        master_pairs = self._get_gene_pairs(fusion_data['master_candidates'])
        all_pairs = self._get_gene_pairs(fusion_data['all_candidates'])
        
        # Get gene networks
        master_networks = self._get_gene_networks(master_pairs)
        all_networks = self._get_gene_networks(all_pairs)
        
        # Add summary section
        self.elements.append(Paragraph("Summary", self.styles.styles["Heading2"]))
        
        # Create summary text
        summary_text = []
        if master_networks:
            summary_text.append(f"Found {len(master_networks)} fusion networks within the target panel:")
            for network in master_networks:
                summary_text.append(f"• {' - '.join(network)}")
        else:
            summary_text.append("No fusion networks identified within the target panel.")
            
        if all_networks and len(all_networks) > len(master_networks):
            summary_text.append(f"\nAdditional {len(all_networks) - len(master_networks)} genome-wide fusion networks identified.")
        
        summary_paragraph = '\n'.join(summary_text)
        self.elements.append(Paragraph(summary_paragraph, self.styles.styles["Normal"]))
        self.elements.append(Spacer(1, 12))
        
        # Add summary to summary section of report
        self.summary_elements.append(Paragraph("Gene Fusions", self.styles.styles["Heading3"]))
        self.summary_elements.append(Paragraph(summary_paragraph, self.styles.styles["Normal"]))
        self.summary_elements.append(Spacer(1, 12))
        
        # Add detailed analysis section if there are any candidates
        if master_pairs or all_pairs:
            self.elements.append(Paragraph("Detailed Analysis", self.styles.styles["Heading2"]))
            self.elements.append(Spacer(1, 12))
            
            # Add target panel fusions
            if master_pairs:
                self.elements.append(Paragraph("Target Panel Fusions", self.styles.styles["Heading3"]))
                table = self._format_fusion_table(fusion_data['master_candidates'], "Target Panel Fusions")
                if table:
                    self.elements.append(table)
                    self.elements.append(Spacer(1, 12))
            
            # Add genome-wide fusions
            if all_pairs and len(all_pairs) > len(master_pairs):
                self.elements.append(Paragraph("Genome-wide Fusions", self.styles.styles["Heading3"]))
                table = self._format_fusion_table(fusion_data['all_candidates'], "Genome-wide Fusions")
                if table:
                    self.elements.append(table)
                    self.elements.append(Spacer(1, 12))
            
            # Add source information
            source_text = "Data sources:\n"
            source_text += f"- Target panel fusions: {fusion_data['master_path']}\n"
            source_text += f"- Genome-wide fusions: {fusion_data['all_path']}\n"
            source_text += "\nNote: Fusion candidates are identified from reads with supplementary alignments and filtered based on mapping quality (>40) and alignment length (>100bp)."
            
            self.elements.append(Paragraph(source_text, self.styles.styles["Normal"]))
        else:
            self.elements.append(Paragraph("No fusion candidates were identified.", self.styles.styles["Normal"]))
