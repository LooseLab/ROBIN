"""
Variant Analysis Section for ROBIN Reports.

This module handles the Pathogenic Variant analysis section of the report.
"""

import os
import logging
import re
import pandas as pd
import numpy as np
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.colors import HexColor
from ..sections.base import ReportSection
from robin import resources

logger = logging.getLogger(__name__)


class VariantResult:
    """Class to store variant analysis results."""
    
    def __init__(self):
        """Initialize variant results."""
        self.snp_data = []
        self.indel_data = []
        self.affected_genes = set()


class VariantAnalysis:
    """Class to handle variant analysis logic."""
    
    def __init__(self):
        """Initialize the variant analyzer."""
        self.gene_info = None
        self._load_resources()
    
    def _load_resources(self):
        """Load required resource files."""
        try:
            gene_info_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "gene_info.txt"
            )
            if os.path.exists(gene_info_file):
                self.gene_info = pd.read_csv(gene_info_file, sep="\t")
                logger.debug("Loaded gene info file with shape: %s", self.gene_info.shape)
        except Exception as e:
            logger.error("Error loading resources: %s", str(e))
            logger.debug("Exception details:", exc_info=True)

    def _is_pathogenic(self, info_str):
        """Check if a variant is pathogenic based on CLNSIG and CLNSIGCONF fields."""
        logger.debug("Checking pathogenicity for variant with info: %s", info_str)
        
        if "CLNSIG=" not in info_str:
            logger.debug("No CLNSIG field found in variant")
            return False
        
        pathogenic_terms = [
            "pathogenic",
            "likely_pathogenic",
            "pathogenic/likely_pathogenic",
            "likely pathogenic",
            "pathogenic/likely pathogenic"
        ]
        
        # First check CLNSIG field
        for field in info_str.split(";"):
            if field.startswith("CLNSIG="):
                clnsig_value = field.split("=")[1].lower()
                logger.debug("Found CLNSIG value: %s", clnsig_value)
                
                # Direct pathogenic classification
                if any(term in clnsig_value for term in pathogenic_terms):
                    logger.debug("Direct pathogenic classification found")
                    return True
                    
                # Handle conflicting classifications
                if "conflicting_classifications_of_pathogenicity" in clnsig_value:
                    logger.debug("Found conflicting classifications, checking CLNSIGCONF")
                    # Look for CLNSIGCONF field
                    for conf_field in info_str.split(";"):
                        if conf_field.startswith("CLNSIGCONF="):
                            conf_value = conf_field.split("=")[1].lower()
                            logger.debug("Found CLNSIGCONF value: %s", conf_value)
                            
                            # Count pathogenic vs benign classifications
                            pathogenic_count = sum(
                                int(count.strip("()")) 
                                for term in pathogenic_terms 
                                for count in conf_value.split("|") 
                                if term in count.lower()
                            )
                            benign_count = sum(
                                int(count.strip("()")) 
                                for count in conf_value.split("|") 
                                if "benign" in count.lower()
                            )
                            logger.debug("Pathogenic count: %d, Benign count: %d", 
                                       pathogenic_count, benign_count)
                            
                            # Return true if there are more pathogenic classifications
                            return pathogenic_count > benign_count
        
        logger.debug("No pathogenic classification found")
        return False

    def process_vcf(self, vcf_file, variant_type, result):
        """Process a VCF file and extract pathological variants."""
        logger.debug("Processing %s VCF file: %s", variant_type, vcf_file)
        try:
            variant_count = 0
            pathogenic_count = 0
            
            with open(vcf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    variant_count += 1
                    if variant_count % 1000 == 0:
                        logger.debug("Processed %d variants", variant_count)
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 8:
                        continue
                    
                    info_str = fields[7]
                    if self._is_pathogenic(info_str):
                        pathogenic_count += 1
                        logger.debug("Found pathogenic variant at %s:%s", fields[0], fields[1])
                        
                        variant_data = {
                            'chromosome': fields[0],
                            'position': int(fields[1]),
                            'reference': fields[3],
                            'alternate': fields[4],
                            'type': variant_type
                        }
                        
                        # Extract additional annotations
                        for field in info_str.split(';'):
                            if field.startswith('GENEINFO='):
                                gene = field.split('=')[1].split(':')[0]
                                variant_data['gene'] = gene
                                result.affected_genes.add(gene)
                                logger.debug("Variant affects gene: %s", gene)
                            elif field.startswith('CLNSIG='):
                                variant_data['significance'] = field.split('=')[1]
                                logger.debug("Clinical significance: %s", variant_data['significance'])
                            elif field.startswith('CLNDN='):
                                variant_data['disease'] = field.split('=')[1]
                                logger.debug("Associated disease: %s", variant_data['disease'])
                        
                        if variant_type == "SNP":
                            result.snp_data.append(variant_data)
                        else:
                            result.indel_data.append(variant_data)
            
            logger.info("Finished processing %s file. Total variants: %d, Pathogenic: %d", 
                       variant_type, variant_count, pathogenic_count)
                            
        except Exception as e:
            logger.error(f"Error processing {variant_type} file {vcf_file}: {str(e)}")
            logger.debug("Exception details:", exc_info=True)


class VariantsSection(ReportSection):
    """Section containing the pathological variants analysis."""

    def __init__(self, report):
        """Initialize the variants section."""
        super().__init__(report)
        self.variant_result = VariantResult()
        self.variant_analyzer = VariantAnalysis()
        self._process_variant_files()

    def setup(self, report):
        """Set up the variants section with the report context."""
        super().setup(report)
        self.variant_result = VariantResult()
        self.variant_analyzer = VariantAnalysis()
        self._process_variant_files()

    def _process_variant_files(self):
        """Process the variant files."""
        logger.debug("Starting variant file processing")
        try:
            # Paths to the annotated VCF files
            snp_vcf = os.path.join(self.report.output, 'clair3', 'snpsift_output.vcf')
            indel_vcf = os.path.join(self.report.output, 'clair3', 'snpsift_indel_output.vcf')
            
            logger.debug("Looking for SNP VCF file at: %s", snp_vcf)
            logger.debug("Looking for indel VCF file at: %s", indel_vcf)
            
            # Process SNPs
            if os.path.exists(snp_vcf):
                logger.debug("Processing SNP VCF file: %s", snp_vcf)
                self.variant_analyzer.process_vcf(snp_vcf, "SNP", self.variant_result)
            else:
                logger.warning("SNP VCF file not found: %s", snp_vcf)
            
            # Process indels
            if os.path.exists(indel_vcf):
                logger.debug("Processing indel VCF file: %s", indel_vcf)
                self.variant_analyzer.process_vcf(indel_vcf, "INDEL", self.variant_result)
            else:
                logger.warning("Indel VCF file not found: %s", indel_vcf)
                
            logger.info("Variant processing complete. Found %d pathogenic SNPs and %d pathogenic indels affecting %d genes",
                       len(self.variant_result.snp_data),
                       len(self.variant_result.indel_data),
                       len(self.variant_result.affected_genes))
                
        except Exception as e:
            logger.error("Error processing variant files: %s", str(e))
            logger.debug("Exception details:", exc_info=True)

    def _format_disease_name(self, disease):
        """Format disease name for better readability."""
        if not disease or disease == 'not_specified':
            return 'Not specified'
        
        # Split on separator characters and clean up
        parts = re.split('[|_]', disease)
        formatted_parts = []
        
        for part in parts:
            # Capitalize first letter of each word
            words = part.split()
            formatted_words = []
            for word in words:
                if word.lower() in ['and', 'of', 'the', 'with', 'or']:
                    formatted_words.append(word.lower())
                else:
                    formatted_words.append(word.capitalize())
            formatted_parts.append(' '.join(formatted_words))
        
        return '\n'.join(formatted_parts)

    def add_content(self):
        """Add the pathogenic variants section to the report."""
        logger.debug("Starting variants section content generation")
        
        # Add page break before detailed section
        self.elements.append(PageBreak())
        
        # Add section title
        self.elements.append(Paragraph("Candidate Pathogenic Variants", self.styles.styles["Heading1"]))
        
        # Add summary section
        self.elements.append(Paragraph("Summary", self.styles.styles["Heading2"]))
        summary_text = (
            f"Found {len(self.variant_result.snp_data)} candidate pathogenic SNPs and "
            f"{len(self.variant_result.indel_data)} candidate pathogenic indels affecting "
            f"{len(self.variant_result.affected_genes)} genes.<br/>"
            f"Genes with candidate pathogenic variants: {', '.join(sorted(self.variant_result.affected_genes))}"
        )
        self.elements.append(Paragraph(summary_text, self.styles.styles["Normal"]))
        self.elements.append(Spacer(1, 12))
        
        # Add summary to summary section
        self.summary_elements.append(Paragraph("Candidate Pathogenic Variants", self.styles.styles["Heading3"]))
        self.summary_elements.append(Paragraph(summary_text, self.styles.styles["Normal"]))
        
        
        if not self.variant_result.snp_data and not self.variant_result.indel_data:
            self.elements.append(Paragraph("No pathogenic variants were identified.", self.styles.styles["Normal"]))
            return
        
        # Add detailed section
        self.elements.append(Paragraph("Detailed Analysis", self.styles.styles["Heading2"]))
        self.elements.append(Spacer(1, 12))
        
        # Create table for pathogenic variants
        if self.variant_result.snp_data or self.variant_result.indel_data:
            # Create header style
            header_style = ParagraphStyle(
                'HeaderStyle',
                parent=self.styles.styles["Normal"],
                fontName='Helvetica-Bold',
                fontSize=8,
                textColor=self.styles.COLORS["primary"],
                alignment=1  # Center alignment
            )
            
            # Create cell style
            cell_style = ParagraphStyle(
                'CellStyle',
                parent=self.styles.styles["Normal"],
                fontSize=8,
                leading=10,  # Line spacing
                spaceBefore=2,
                spaceAfter=2
            )
            
            # Format headers as paragraphs
            table_data = [[
                Paragraph('Type', header_style),
                Paragraph('Chr', header_style),
                Paragraph('Position', header_style),
                Paragraph('Gene', header_style),
                Paragraph('Change', header_style),
                Paragraph('Disease', header_style)
            ]]
            
            # Add variants to table
            for variant_list in [self.variant_result.snp_data, self.variant_result.indel_data]:
                for variant in variant_list:
                    change = f"{variant['reference']}>{variant['alternate']}"
                    disease = self._format_disease_name(variant.get('disease', 'Not specified'))
                    
                    row = [
                        Paragraph(variant['type'], cell_style),
                        Paragraph(variant['chromosome'], cell_style),
                        Paragraph(str(variant['position']), cell_style),
                        Paragraph(variant.get('gene', 'Unknown'), cell_style),
                        Paragraph(change, cell_style),
                        Paragraph(disease, cell_style)
                    ]
                    table_data.append(row)
            
            # Adjust column widths (total should be around 7 inches for A4 paper)
            col_widths = [inch * x for x in [0.5, 0.4, 0.8, 0.8, 0.7, 3.0]]
            
            # Create and style the table
            table = Table(table_data, colWidths=col_widths, repeatRows=1)
            table.setStyle(TableStyle([
                # Inherit modern table style
                *self.MODERN_TABLE_STYLE._cmds,
                
                # Preserve specific alignments
                ('ALIGN', (0, 0), (0, -1), 'CENTER'),  # Type column centered
                ('ALIGN', (1, 0), (1, -1), 'CENTER'),  # Chr column centered
                ('ALIGN', (2, 0), (2, -1), 'RIGHT'),   # Position column right-aligned
                ('ALIGN', (3, 0), (3, -1), 'LEFT'),    # Gene column left-aligned
                ('ALIGN', (4, 0), (4, -1), 'CENTER'),  # Change column centered
                ('ALIGN', (5, 0), (5, -1), 'LEFT'),    # Disease column left-aligned
            ]))
            
            self.elements.append(table)
            self.elements.append(Spacer(1, 6))
            
            # Add source file information
            source_style = ParagraphStyle(
                'Source',
                parent=self.styles.styles["Normal"],
                fontSize=8,
                textColor=HexColor("#6B7280"),
                leading=10
            )
            
            # Get the VCF file paths
            snp_vcf = os.path.join(self.report.output, 'clair3', 'snpsift_output.vcf')
            indel_vcf = os.path.join(self.report.output, 'clair3', 'snpsift_indel_output.vcf')
            
            # Format the source text
            source_text = "Data sources:"
            if os.path.exists(snp_vcf):
                source_text += f"\nSNPs: {snp_vcf}"
            if os.path.exists(indel_vcf):
                source_text += f"\nIndels: {indel_vcf}"
            
            self.elements.append(Paragraph(source_text, source_style))
            self.elements.append(Spacer(1, 12))
            
            # Add a note about the variants
            note_style = ParagraphStyle(
                'Note',
                parent=self.styles.styles["Normal"],
                fontSize=8,
                textColor=HexColor("#4B5563"),
                leading=10
            )
            note_text = (
                "Note: Variants are classified as pathogenic based on ClinVar annotations. "
                "Disease associations are derived from ClinVar's CLNDN field where available."
            )
            self.elements.append(Paragraph(note_text, note_style))
        
        self.elements.append(PageBreak()) 