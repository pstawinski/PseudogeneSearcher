package pl.genebeam.pseudogenes.helpers;

import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.google.common.base.Splitter;
import com.google.common.collect.Range;
import com.google.common.collect.TreeRangeMap;

public class Transcript {
	private String name;
	private String gene;
	private String chr;
	private int txStart;
	private int txEnd;
	private int cdsStart;
	private int cdsEnd;
	@JsonIgnore
	private List<Integer> exonStarts;
	@JsonIgnore
	private List<Integer> intronStarts;
	private boolean strandPlus;

	/**
	 * Range -> exone index
	 */
	@JsonIgnore
	transient private TreeRangeMap<Integer, Integer> exoneRangeMap, introneRangeMap;
	@JsonIgnore
	transient private TreeSet<Integer> junctionSites;
	@JsonIgnore
	transient private TreeSet<Integer> exonStartPointSet, intronStartPointSet;

	private final static Splitter tabSplitter = Splitter.on('\t');
	private final static Splitter commaSplitter = Splitter.on(',').trimResults().omitEmptyStrings();

	public Transcript(String name, String gene, String chr, int txStart, int txEnd, int cdsStart, int cdsEnd,
			List<Integer> exonStarts, List<Integer> intronStarts, boolean strandPlus) {
		super();
		this.name = name;
		this.gene = gene;
		this.chr = chr;
		this.txStart = txStart;
		this.txEnd = txEnd;
		this.cdsStart = cdsStart;
		this.cdsEnd = cdsEnd;
		this.exonStarts = exonStarts;
		this.intronStarts = intronStarts;
		this.strandPlus = strandPlus;
	}

	public Transcript(String line) {
		List<String> lineSplitted = tabSplitter.splitToList(line);
		this.name = lineSplitted.get(0);
		this.chr = lineSplitted.get(1);
		this.strandPlus = "+".equals(lineSplitted.get(2));
		this.txStart = Integer.valueOf(lineSplitted.get(3));
		this.txEnd = Integer.valueOf(lineSplitted.get(4));
		this.cdsStart = Integer.valueOf(lineSplitted.get(5));
		this.cdsEnd = Integer.valueOf(lineSplitted.get(6));
		// this.exonCount = Integer.valueOf(lineSplitted.get(7));
		this.exonStarts = commaSplitter.splitToList(lineSplitted.get(8)).stream().map(Integer::valueOf).map(x -> x + 1)
				.collect(Collectors.toList());
		this.intronStarts = commaSplitter.splitToList(lineSplitted.get(9)).stream().map(Integer::valueOf)
				.map(x -> x + 1).collect(Collectors.toList());

		// this.score = lineSplitted.get(10);
		this.gene = lineSplitted.get(11);
		// this.cdsStartStat = lineSplitted.get(12);
		// this.cdsEndStat = lineSplitted.get(13);
	}

	@JsonIgnore
	public int getCdsEnd() {
		return cdsEnd;
	}

	@JsonIgnore
	public int getCdsStart() {
		return cdsStart;
	}

	@JsonIgnore
	public List<Integer> getIntronStarts() {
		return intronStarts;
	}

	@JsonIgnore
	public List<Integer> getExonStarts() {
		return exonStarts;
	}

	public String getGene() {
		return gene;
	}

	public String getName() {
		return name;
	}

	@JsonIgnore
	public int getTxEnd() {
		return txEnd;
	}

	@JsonIgnore
	public int getTxStart() {
		return txStart;
	}

	public String getChr() {
		return chr;
	}

	@JsonIgnore
	public TreeRangeMap<Integer, Integer> getExoneMap() {
		if (exoneRangeMap == null) {
			exoneRangeMap = TreeRangeMap.create();
			for (int i = 0; i < exonStarts.size(); i++) {
				exoneRangeMap.put(Range.closedOpen(exonStarts.get(i), intronStarts.get(i)), i);
			}
		}
		return exoneRangeMap;
	}

	@JsonIgnore
	public TreeRangeMap<Integer, Integer> getIntroneMap() {

		if (introneRangeMap == null) {
			introneRangeMap = TreeRangeMap.create();
			for (int i = 0; i < exonStarts.size() - 1; i++) {
				introneRangeMap.put(Range.closedOpen(intronStarts.get(i), exonStarts.get(i + 1)), i);
			}
		}

		return introneRangeMap;
	}

	public Integer getExonIndexAtPosition(int position) {
		return getExoneMap().get(position);
	}

	public Integer getIntronIndexAtPosition(int position) {
		return getIntroneMap().get(position);
	}

	public Collection<Integer> getExonesIn(Range<Integer> range) {
		return getExoneMap().subRangeMap(range).asMapOfRanges().values();
	}

	/**
	 * get intrones overlapping in this range
	 * 
	 * @param range
	 * @return
	 */
	public Collection<Integer> getIntronesOverlapping(Range<Integer> range) {
		return getIntroneMap().subRangeMap(range).asMapOfRanges().values();
	}

	public Collection<Integer> getIntronesIn(Range<Integer> range) {
		return getIntroneMap().asMapOfRanges().entrySet().stream().filter((x) -> range.encloses(x.getKey()))
				.map(x -> x.getValue()).collect(Collectors.toList());
	}

	public Integer getIntronesSizesIn(Range<Integer> range) {
		return getIntroneMap().subRangeMap(range).asMapOfRanges().keySet().stream()
				.mapToInt(x -> x.upperEndpoint() - x.lowerEndpoint()).sum();
	}

	public Range<Integer> getIntronWithIndex(int index) {
		return Range.closedOpen(intronStarts.get(index), exonStarts.get(index + 1));
	}

	public Range<Integer> getExonWithIndex(int index) {
		return Range.closedOpen(exonStarts.get(index), intronStarts.get(index));
	}

	public TreeSet<Integer> getJunctionSites() {
		if (junctionSites == null) {
			junctionSites = new TreeSet<>();
			junctionSites.addAll(exonStarts);
			junctionSites.addAll(intronStarts);
		}
		return junctionSites;
	}

	public Integer getNearestJunctionSite(Integer position) {
		TreeSet<Integer> js = getJunctionSites();
		Integer ceiling = Optional.ofNullable(js.ceiling(position)).orElse(Integer.MAX_VALUE);
		Integer floor = Optional.ofNullable(js.floor(position)).orElse(0);

		return Math.min(ceiling - position, position - floor);

	}

	public SortedSet<Integer> getJunctionsIn(Range<Integer> range) {
		TreeSet<Integer> js = getJunctionSites();
		return js.subSet(range.lowerEndpoint(), range.upperEndpoint());
	}

	@JsonIgnore
	public SortedSet<Integer> getExonStartPointsSet() {
		if (exonStartPointSet == null) {
			exonStartPointSet = new TreeSet<>(exonStarts);
		}
		return exonStartPointSet;
	}

	@JsonIgnore
	public SortedSet<Integer> getIntronStartPointsSet() {
		if (intronStartPointSet == null) {
			intronStartPointSet = new TreeSet<>(intronStarts);
		}
		return intronStartPointSet;
	}

	public Boolean isExonStartPoint(Integer position) {
		return getExonStartPointsSet().contains(position);
	}

	public Boolean isIntronStartPoint(Integer position) {
		return getIntronStartPointsSet().contains(position);
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Transcript other = (Transcript) obj;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		return true;
	}

	@JsonIgnore
	public int getIntronesNumber() {
		return intronStarts.size() - 1;
	}

	@Override
	public String toString() {
		return getGene() + ":" + getName();
	}
}
