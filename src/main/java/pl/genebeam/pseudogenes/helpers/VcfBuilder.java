package pl.genebeam.pseudogenes.helpers;

import java.io.Closeable;
import java.io.OutputStream;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import pl.genebeam.pseudogenes.model.IdentifiedPseudogene;

public class VcfBuilder implements Closeable {
	private final VariantContextWriter vcf;
	private final ReferenceSequenceFile reference;

	public VcfBuilder(OutputStream os, ReferenceSequenceFile reference, String sampleName) {
		this.vcf = new VariantContextWriterBuilder().setOutputStream(os)
				.setOptions(EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER)).build();
		this.reference = reference;
		buildHeader(sampleName);
	}

	private void buildHeader(String sampleName) {
		Set<VCFHeaderLine> metadata = new HashSet<>();

		metadata.add(
				new VCFFormatHeaderLine("ITRSP", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "Number of ... "));
		metadata.add(
				new VCFInfoHeaderLine("AvaStr", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Id blabla"));

		// VCFHeader header = new Vcfh;
		VCFHeader header = new VCFHeader(metadata, Collections.singleton(sampleName));
		// vcf.setHeader(header);
		vcf.writeHeader(header);
	}

	public void add(IdentifiedPseudogene pseudogene, String sampleName) {
		if (pseudogene.getFracIntronesCoveredClip() < 0.1 && pseudogene.getFracIntronesCoveredIS() < 0.1
				&& pseudogene.getFracIntronesCoveredCombined() < 0.1) {
			// not enough data to say anything, do not save
		} else {
			int start = pseudogene.getMotherTranscript().getCdsStart();
			int end = pseudogene.getMotherTranscript().getCdsEnd();

			String refBase = reference.getSubsequenceAt(pseudogene.getContig(), start, start).getBaseString();

			VariantContextBuilder vcBuilder = new VariantContextBuilder().chr(pseudogene.getContig())
					.start(pseudogene.getStart()).stop(end).attribute("SVTYPE", "PSDGN").attribute("END", end)
					.attribute("GENE", pseudogene.getMotherTranscript().getGene())
					.attribute("TRANSCRIPT", pseudogene.getMotherTranscript().getName()).alleles(refBase, "<PDG>");

			if (pseudogene.isReferencePseudogene()) {
				vcBuilder.attribute("REFSDGN", true);
			}
			if (pseudogene.isReferencePseudogene()) {
				vcBuilder.attribute("ABYSDGN",
						pseudogene.getAbysovPseudogene().stream().collect(Collectors.joining(",")));
			}

			if (pseudogene.isPass()) {
				vcBuilder.passFilters();
			} else {
				vcBuilder.filters(pseudogene.getFilters());
			}

			List<IntronData> intronDetails = pseudogene.getIntronDetails();
			vcBuilder.attribute("INTRONES", pseudogene.getMotherTranscript().getIntronesNumber());

			GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName,
					Collections.singletonList(Allele.create("<PDG>".getBytes(), false)));
			genotypeBuilder.attribute("CL", pseudogene.getFracIntronesCoveredClip());
			genotypeBuilder.attribute("IS", pseudogene.getFracIntronesCoveredIS());
			genotypeBuilder.attribute("CO", pseudogene.getFracIntronesCoveredCombined());

			vcBuilder.genotypes(genotypeBuilder.make());

			vcf.add(vcBuilder.make());
		}
	}

	public void close() {
		vcf.close();
	}
}
