package pl.genebeam.pseudogenes.model;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.annotation.Nullable;

import com.fasterxml.jackson.annotation.JsonIgnore;

import pl.genebeam.pseudogenes.helpers.IntronData;
import pl.genebeam.pseudogenes.helpers.Transcript;

public class IdentifiedPseudogene {
	private String contig;
	/**
	 * position of mother gene
	 */
	private int start;
	@Nullable
	private Set<String> filters = new HashSet<>();
	private boolean pass = false;

	private Transcript motherTranscript;
	@Nullable

	private boolean referencePseudogene;
	@JsonIgnore
	private Collection<String> abysovPseudogene;

	@JsonIgnore
	private Collection<String> disconcordances;
	private List<IntronData> intronDetails;
	private boolean known;
	private double fracIntronesCoveredIS;
	private double fracIntronesCoveredClip;
	@JsonIgnore
	private double fracIntronesCoveredCombined;
	private String sampleName;

	private List<PositionWithEvidence> possiblePositions;

	public String getContig() {
		return contig;
	}

	public void setContig(String contig) {
		this.contig = contig;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public Set<String> getFilters() {
		return filters;
	}

	public void setFilters(Set<String> filters) {
		this.filters = filters;
	}

	public boolean isPass() {
		return pass;
	}

	public void setPass(boolean pass) {
		this.pass = pass;
	}

	public Transcript getMotherTranscript() {
		return motherTranscript;
	}

	public void setMotherTranscript(Transcript motherTranscript) {
		this.motherTranscript = motherTranscript;
	}

	public void setReferencePseudogene(boolean b) {
		this.referencePseudogene = b;

	}

	public void setAbysovPseudogene(Collection<String> collection) {
		this.abysovPseudogene = collection;
	}

	public boolean isReferencePseudogene() {
		return referencePseudogene;
	}

	public void setFracIntronesCoveredCombined(double fracIntronesCoveredCombined) {
		this.fracIntronesCoveredCombined = fracIntronesCoveredCombined;

	}

	public void setFracIntronesCoveredClip(double fracIntronesCoveredClip) {
		this.fracIntronesCoveredClip = fracIntronesCoveredClip;

	}

	public void setFracIntronesCoveredIS(double fracIntronesCoveredIS) {
		this.fracIntronesCoveredIS = fracIntronesCoveredIS;

	}

	public void setKnown(boolean known) {
		this.known = known;

	}

	public void setIntronDetails(List<IntronData> intronDetails) {
		this.intronDetails = intronDetails;

	}

	public void setDisconcordances(Collection<String> disconcordances) {
		this.disconcordances = disconcordances;

	}

	public Collection<String> getAbysovPseudogene() {
		return abysovPseudogene;
	}

	public Collection<String> getDisconcordances() {
		return disconcordances;
	}

	public List<IntronData> getIntronDetails() {
		return intronDetails;
	}

	public boolean isKnown() {
		return known;
	}

	public double getFracIntronesCoveredIS() {
		return fracIntronesCoveredIS;
	}

	public double getFracIntronesCoveredClip() {
		return fracIntronesCoveredClip;
	}

	public double getFracIntronesCoveredCombined() {
		return fracIntronesCoveredCombined;
	}

	public List<PositionWithEvidence> getPossiblePositions() {
		return possiblePositions;
	}

	public void setPossiblePositions(List<PositionWithEvidence> possiblePositions) {
		this.possiblePositions = possiblePositions;
	}

	public String getSampleName() {
		return sampleName;
	}

	public void setSampleName(String sampleName) {
		this.sampleName = sampleName;
	}

}
