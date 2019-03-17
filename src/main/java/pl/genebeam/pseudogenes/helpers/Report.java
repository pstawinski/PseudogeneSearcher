package pl.genebeam.pseudogenes.helpers;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.collections4.CollectionUtils;

import com.fasterxml.jackson.core.JsonFactory;
import com.fasterxml.jackson.core.JsonGenerationException;
import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.collect.ComparisonChain;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import pl.genebeam.pseudogenes.model.AlignmentResult;
import pl.genebeam.pseudogenes.model.IdentifiedPseudogene;
import pl.genebeam.pseudogenes.model.PositionWithEvidence;
import pl.genebeam.pseudogenes.model.ReadSummary;
import pl.genebeam.utils.TxNameToGeneName;

public class Report {
	private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(Report.class);

	private final static NumberFormat formatter = new DecimalFormat("#0.00");
	private Map<Transcript, TranscriptStats> transcriptStats = new HashMap<>();
	private Map<Transcript, Collection<String>> disconcordances = new HashMap<>();

	private Map<Transcript, Collection<String>> disconcordancesAbysov = new HashMap<>();
	private Map<Transcript, Boolean> abysovGene = new HashMap<>();

	private final String TX_NAME_DEBUG_POINT = "NR_026803";
	private Transcript debugTranscript;

	private TxNameToGeneName txNameToGeneName;

	private Set<Transcript> referencePseudogenes = new HashSet<>();

	private Map<Transcript, List<PositionWithEvidence>> possiblePositions = new HashMap<>();

	public Report(TxNameToGeneName txNameToGeneName) {
		this.txNameToGeneName = txNameToGeneName;
	}

	public String toString() {
		List<Entry<Transcript, TranscriptStats>> sorted = transcriptStats.entrySet().stream()
				.sorted((a, b) -> ComparisonChain.start()
						.compare(a.getValue().getFracIntronesCoveredCombined(),
								b.getValue().getFracIntronesCoveredCombined())
						.compare(a.getValue().getNumberOfIntrons(), b.getValue().getNumberOfIntrons()).result())
				.collect(Collectors.toList());

		StringBuilder sb = new StringBuilder();

		// remove duplicated, leave only the most important one
		Collections.reverse(sorted);
		Set<String> geneNameUniqueness = new HashSet<>();
		for (Iterator<Entry<Transcript, TranscriptStats>> iterator = sorted.iterator(); iterator.hasNext();) {
			Entry<Transcript, TranscriptStats> entry = iterator.next();
			if (geneNameUniqueness.add(entry.getKey().getGene())) {
				// it's ok, we see this gene for the first time
			} else {
				iterator.remove();
			}
		}

		Collections.reverse(sorted);
		for (Entry<Transcript, TranscriptStats> tsStats : sorted) {

			boolean known = false;
			if (CollectionUtils.isNotEmpty(disconcordances.get(tsStats.getKey()))
					|| CollectionUtils.isNotEmpty(disconcordancesAbysov.get(tsStats.getKey()))
					|| disconcordancesAbysov.get(tsStats.getKey()) != null
					|| referencePseudogenes.contains(tsStats.getKey()))
				known = true;

			sb.append(tsStats.getKey().getGene() + ":" + tsStats.getKey().getName());
			sb.append("\t");
			sb.append(known ? "KNOWN" : "");
			sb.append("\t" + tsStats.getKey().getChr() + ":" + tsStats.getKey().getTxStart() + "-"
					+ tsStats.getKey().getTxEnd() + "\t");
			sb.append("\t" + formatter.format(tsStats.getValue().getFracIntronesCoveredCombined()) + "\t");
			sb.append("\t" + formatter.format(tsStats.getValue().getFracIntronesCoveredClip()) + "\t");
			sb.append("\t" + formatter.format(tsStats.getValue().getFracIntronesCoveredIS()) + "\t");

			tsStats.getValue().supportingCount().stream().forEach(intron -> {
				sb.append(" (");
				sb.append(intron.getIndex() + ":");
				sb.append(intron.getReadsSupportingRemovalByAlignment() + ",");
				sb.append(intron.getReadsSupportingRemovalByIS() + ")");
			});

			sb.append("\tdisconcordances:" + disconcordances.get(tsStats.getKey()));
			sb.append("\tabysov:" + disconcordancesAbysov.get(tsStats.getKey()));
			sb.append("\tabysovGene:" + (disconcordancesAbysov.get(tsStats.getKey())) != null ? "KNOWN" : "");
			sb.append("\tpseudoGene:" + (referencePseudogenes.contains(tsStats.getKey()) ? "KNOWN" : ""));

			sb.append("\n");
		}

		if (TX_NAME_DEBUG_POINT != null) {
			TranscriptStats debugTranscriptSTats = getTranscriptStats(debugTranscript, false);
			if (debugTranscriptSTats != null)
				System.out.println(debugTranscriptSTats);
		}

		return sb.toString();

	}

	public void toVcf(ReferenceSequenceFile reference, OutputStream os, OutputStream jsonOs, String sampleName)
			throws JsonGenerationException, JsonMappingException, IOException {
		VcfBuilder vcfBuilder = new VcfBuilder(os, reference, sampleName);

		List<Entry<Transcript, TranscriptStats>> sorted = transcriptStats.entrySet().stream()
				.sorted((a, b) -> ComparisonChain.start()
						.compare(a.getValue().getFracIntronesCoveredCombined(),
								b.getValue().getFracIntronesCoveredCombined())
						.compare(a.getValue().getNumberOfIntrons(), b.getValue().getNumberOfIntrons()).result())
				.collect(Collectors.toList());

		StringBuilder sb = new StringBuilder();

		// remove duplicated, leave only the most important one
		Collections.reverse(sorted);
		Set<String> geneNameUniqueness = new HashSet<>();
		for (Iterator<Entry<Transcript, TranscriptStats>> iterator = sorted.iterator(); iterator.hasNext();) {
			Entry<Transcript, TranscriptStats> entry = iterator.next();
			if (geneNameUniqueness.add(entry.getKey().getGene())) {
				// it's ok, we see this gene for the first time
			} else {
				iterator.remove();
			}
		}

		JsonFactory jsonFactory = new JsonFactory();
		jsonFactory.configure(JsonGenerator.Feature.AUTO_CLOSE_TARGET, false);
		ObjectMapper objectMapper = new ObjectMapper(jsonFactory);

		Collections.reverse(sorted);
		for (Entry<Transcript, TranscriptStats> tsStats : sorted) {
			IdentifiedPseudogene ip = new IdentifiedPseudogene();

			ip.setMotherTranscript(tsStats.getKey());
			ip.setContig(tsStats.getKey().getChr());
			ip.setStart(tsStats.getKey().getTxStart());
			ip.setSampleName(sampleName);

			boolean known = false;
			if (CollectionUtils.isNotEmpty(disconcordances.get(tsStats.getKey()))
					|| CollectionUtils.isNotEmpty(disconcordancesAbysov.get(tsStats.getKey()))
					|| disconcordancesAbysov.get(tsStats.getKey()) != null
					|| referencePseudogenes.contains(tsStats.getKey())) {
				known = true;
			}
			ip.setKnown(known);
			ip.setFracIntronesCoveredCombined(tsStats.getValue().getFracIntronesCoveredCombined());
			ip.setFracIntronesCoveredClip(tsStats.getValue().getFracIntronesCoveredClip());
			ip.setFracIntronesCoveredIS(tsStats.getValue().getFracIntronesCoveredIS());

			// sb.append("\t" +
			// formatter.format(tsStats.getValue().getFracIntronesCoveredCombined())
			// + "\t");
			// sb.append("\t" +
			// formatter.format(tsStats.getValue().getFracIntronesCoveredClip())
			// + "\t");
			// sb.append("\t" +
			// formatter.format(tsStats.getValue().getFracIntronesCoveredIS()) +
			// "\t");

			ip.setIntronDetails(tsStats.getValue().supportingCount());

			// tsStats.getValue().supportingCount().stream().forEach(intron -> {
			// sb.append(" (");
			// sb.append(intron.getIndex() + ":");
			// sb.append(intron.getReadsSupportingRemovalByAlignment() + ",");
			// sb.append(intron.getReadsSupportingRemovalByIS() + ")");
			// });

			// sb.append("\tdisconcordances:" +
			// disconcordances.get(tsStats.getKey()));
			ip.setDisconcordances(disconcordances.get(tsStats.getKey()));
			// sb.append("\tabysov:" +
			// disconcordancesAbysov.get(tsStats.getKey()));
			// sb.append("\tabysovGene:" +
			// (disconcordancesAbysov.get(tsStats.getKey())) != null ? "KNOWN" :
			// "");
			ip.setAbysovPseudogene(disconcordancesAbysov.get(tsStats.getKey()));
			// sb.append("\tpseudoGene:" +
			// (referencePseudogenes.contains(tsStats.getKey()) ? "KNOWN" :
			// ""));
			ip.setReferencePseudogene(referencePseudogenes.contains(tsStats.getKey()));
			sb.append("\n");
			ip.setPossiblePositions(possiblePositions.get(ip.getMotherTranscript()));

			if (ip.getFracIntronesCoveredClip() > 0.5 && ip.getFracIntronesCoveredIS() > 0.5
					&& ip.getFracIntronesCoveredCombined() > 0.5 && ip.getIntronDetails().size() > 1) {
				ip.setPass(true);
			} else {
				ip.setPass(false);

				Set<String> filters = new HashSet<>();
				if (ip.getFracIntronesCoveredClip() <= 0.5) {
					filters.add("LOW_CLIP");
				}
				if (ip.getFracIntronesCoveredIS() <= 0.5) {
					filters.add("LOW_IS");
				}
				if (ip.getFracIntronesCoveredCombined() <= 0.5) {
					filters.add("LOW_COMB");
				}
				if (ip.getIntronDetails().size() <= 1) {
					filters.add("LOW_INTRON_COUNT");
				}

				ip.setFilters(filters);
			}

			vcfBuilder.add(ip, sampleName);
			objectMapper.writeValue(jsonOs, ip);
			jsonOs.write("\n".getBytes(StandardCharsets.UTF_8));
		}

		if (TX_NAME_DEBUG_POINT != null) {
			TranscriptStats debugTranscriptSTats = getTranscriptStats(debugTranscript, false);
			if (debugTranscriptSTats != null)
				System.out.println(debugTranscriptSTats);
		}

		jsonOs.close();
		vcfBuilder.close();

	}

	public void addByInsertSize(Transcript transcript, Integer txIdx, SAMRecord read) {
		if (TX_NAME_DEBUG_POINT.equals(transcript.getName())) {
			log.warn("Here");
			debugTranscript = transcript;
		}

		TranscriptStats ts = getTranscriptStats(transcript, true);
		ts.addReadWithSmallInsertSize(txIdx, read);

	}

	public void addByAlign(Transcript transcript, Integer txIdx, SAMRecord read, AlignmentResult ar) {
		if (TX_NAME_DEBUG_POINT.equals(transcript.getName())) {
			log.warn("Here");
		}

		TranscriptStats ts = getTranscriptStats(transcript, true);
		ts.addReadWithSoftClip(txIdx, read, ar);

	}

	private TranscriptStats getTranscriptStats(Transcript t, boolean create) {
		TranscriptStats ts = transcriptStats.get(t);
		if (ts == null && create) {
			synchronized (transcriptStats) {
				ts = transcriptStats.get(t);
				if (ts == null) {
					ts = new TranscriptStats(t);
					transcriptStats.put(t, ts);
				}
			}
		}
		return ts;
	}

	public void addDisconcordances(Transcript transcript, Set<String> hitsPseudogenesOf) {
		disconcordances.put(transcript, hitsPseudogenesOf);

	}

	public void addAbysovDisconcordances(Transcript transcript, Set<String> hitsAbyzovOf) {
		disconcordancesAbysov.put(transcript, hitsAbyzovOf);

	}

	public void addAbysovGene(Transcript transcript, boolean abysov) {
		abysovGene.put(transcript, abysov);

	}

	public void addReferencePseudogene(Transcript transcript) {
		referencePseudogenes.add(transcript);

	}

	public void addPossiblePositions(Transcript transcript, List<PositionWithEvidence> possiblePositions) {
		this.possiblePositions.put(transcript, possiblePositions);

	}

	private HashMap<String, ReadSummary> readsMap = null;

	public synchronized void populateReadData(SAMRecord read) {
		if (readsMap == null) {
			readsMap = new HashMap<>();
			possiblePositions.values().stream().flatMap(v -> v.stream()).flatMap(v -> v.getReads().stream())
					.forEach(r -> readsMap.put(r.getName(), r));
		}
		if (readsMap.containsKey(read.getReadName())) {
			readsMap.get(read.getReadName()).setBases(read.getReadString());
		}

	}

}
