package pl.genebeam.pseudogenes.helpers;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;

import gnu.trove.set.hash.TIntHashSet;
import htsjdk.samtools.SAMRecord;
import pl.genebeam.pseudogenes.model.AlignmentResult;

public class TranscriptStats {
	public TranscriptStats(Transcript t) {
		this.transcript = t;
	}

	private final Transcript transcript;
	private Multimap<Integer, SAMRecord> readsByIntronSize = Multimaps.synchronizedSetMultimap(HashMultimap.create());
	private Multimap<Integer, AlignmentResult> readsBySoftClipped = Multimaps
			.synchronizedSetMultimap(HashMultimap.create());

	public void addReadWithSmallInsertSize(Integer txIdx, SAMRecord read) {
		readsByIntronSize.put(txIdx, read);
	}

	public void addReadWithSoftClip(Integer txIdx, SAMRecord read, AlignmentResult ar) {
		readsBySoftClipped.put(txIdx, ar);
	}

	public double getFracIntronesCoveredIS() {
		double size = readsByIntronSize.keySet().size();
		return size / (double) transcript.getIntronesNumber();
	}

	public double getFracIntronesCoveredClip() {
		double size = readsBySoftClipped.keySet().size();
		return size / (double) transcript.getIntronesNumber();
	}

	public double getFracIntronesCoveredCombined() {
		TIntHashSet set = new TIntHashSet();
		set.addAll(readsByIntronSize.keySet());
		set.addAll(readsBySoftClipped.keySet());
		// return ((double) set.size()) / (double) getNumberOfIntrons();
		double doubleCoveredFraction = (double) set.size() / (double) transcript.getIntronesNumber();
		return (doubleCoveredFraction + 2 * getFracIntronesCoveredClip() + getFracIntronesCoveredIS()) / 4.;

	}

	public List<IntronData> supportingCount() {
		List<IntronData> list = IntStream.range(0, transcript.getIntronesNumber())
				.filter((x) -> (!readsByIntronSize.get(x).isEmpty()) || (!readsBySoftClipped.get(x).isEmpty()))
				.mapToObj((x) -> new IntronData(x, readsByIntronSize.get(x).size(), readsBySoftClipped.get(x).size(),
						readsBySoftClipped.get(x)))
				.collect(Collectors.toList());
		return list;
	}

	public int getNumberOfIntrons() {
		try {
			return transcript.getIntronesNumber();
		} catch (Exception e) {
			throw e;
		}
	}

}
