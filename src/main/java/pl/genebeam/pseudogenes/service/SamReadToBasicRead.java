package pl.genebeam.pseudogenes.service;

import htsjdk.samtools.SAMRecord;
import pl.genebeam.pseudogenes.model.ReadBasicData;

public class SamReadToBasicRead {
	public static ReadBasicData convert(SAMRecord read) {
		ReadBasicData rbd = new ReadBasicData();
		rbd.setName(read.getReadName());
		rbd.setPosition(read.getAlignmentStart());
		rbd.setContig(read.getContig().intern());		
		return rbd;
	}
}
