package pl.genebeam.pseudogenes.helpers;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.ComparisonChain;

import htsjdk.samtools.SAMRecord;

public class Report {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(Report.class);

    private final static NumberFormat formatter = new DecimalFormat("#0.00");
    private Map<Transcript, TranscriptStats> transcriptStats = new HashMap<>();

    private final String TX_NAME_DEBUG_POINT = "NR_026803";
    private Transcript debugTranscript;

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

            sb.append(tsStats.getKey().getGene() + ":" + tsStats.getKey().getName());
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
            sb.append("\n");
        }

        if (TX_NAME_DEBUG_POINT != null) {
            TranscriptStats debugTranscriptSTats = getTranscriptStats(debugTranscript, false);
            if (debugTranscriptSTats != null)
                System.out.println(debugTranscriptSTats);
        }

        return sb.toString();

    }

    public void addByInsertSize(Transcript transcript, Integer txIdx, SAMRecord read) {
        if (TX_NAME_DEBUG_POINT.equals(transcript.getName())) {
            log.warn("Here");
            debugTranscript = transcript;
        }

        TranscriptStats ts = getTranscriptStats(transcript, true);
        ts.addReadWithSmallInsertSize(txIdx, read);

    }

    public void addByAlign(Transcript transcript, Integer txIdx, SAMRecord read) {
        if (TX_NAME_DEBUG_POINT.equals(transcript.getName())) {
            log.warn("Here");
        }

        TranscriptStats ts = getTranscriptStats(transcript, true);
        ts.addReadWithSoftClip(txIdx, read);

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

}
