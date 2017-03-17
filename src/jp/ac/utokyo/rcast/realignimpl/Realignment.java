package jp.ac.utokyo.rcast.realignimpl;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.StringUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.StringUtil;
import jp.ac.utokyo.rcast.karkinos.exec.DataSet;
import jp.ac.utokyo.rcast.karkinos.utils.OptionComparator;
import jp.ac.utokyo.rcast.karkinos.utils.ReadWriteBase;
import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;
import srma.SRMA_mod;
import srma.TwobitReferenceSequence;

public class Realignment extends ReadWriteBase {

	public static void main(String[] arg) {

		//
		BasicParser parcer = new BasicParser();
		List<Option> optionList = getOptionListForKarkinos();
		Options opts = new Options();
		for (Option opt : optionList) {
			opts.addOption(opt);
		}

		CommandLine cl = null;
		try {
			cl = parcer.parse(opts, arg);
		} catch (ParseException e1) {
			System.out.println(e1.getMessage());
			HelpFormatter help = new HelpFormatter();
			help.setOptionComparator(new OptionComparator(optionList));
			help.printHelp("dualRealign.jar", opts, true);
			return;
		}

		Realignment ra = new Realignment();
		try {
			ra.exec(cl);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private void exec(CommandLine cl) throws Exception {

		String normalbamf = cl.getOptionValue("n");
		String tumorbamf = cl.getOptionValue("t");
		String twobitref = cl.getOptionValue("r");

		String targetRegion = null;
		if (cl.hasOption("ct")) {
			targetRegion = cl.getOptionValue("ct");
		}
		String outdir = cl.getOptionValue("o");
		if (!outdir.endsWith("/")) {
			outdir = outdir + "/";
			File f = new File(outdir);
			if (!f.exists()) {
				boolean suc = false;
				try {
					suc = f.mkdirs();
				} catch (Exception ex) {
					System.out.println("could not make directory " + outdir);
					return;
				}
				if (suc == false) {
					System.out.println("could not make directory " + outdir);
					return;
				}
			}
		}

		int numthread = 1;
		if (cl.hasOption("nt")) {
			numthread = Integer.parseInt(cl.getOptionValue("nt"));
		}

		try {
			realign(normalbamf, tumorbamf, twobitref, targetRegion, outdir, numthread);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private void realign(String normalbamf, String tumorbamf, String twobitref, String targetRegion, String outdir,
			int numthread) throws Exception {

		// load target
		TwoBitGenomeReader tgr = new TwoBitGenomeReader(new File(twobitref));
		tgr.setCheckRepeat(false);
		DataSet dataset = null;
		if (targetRegion != null) {
			dataset = new DataSet(tgr.readIndex());
			dataset.loadTargetBed(targetRegion, tgr);
		}
		//

		SAMFileReader forheader = getReader(normalbamf);
		String outnameN = outdir + "/" + new File(normalbamf).getName() + "_realign.bam";
		String outnameT = outdir + "/" + new File(tumorbamf).getName() + "_realign.bam";

		SAMFileWriter normalbamw = getPreSortWriter(forheader.getFileHeader(), outnameN);
		SAMFileWriter tumorbamw = getPreSortWriter(forheader.getFileHeader(), outnameT);

		List<SAMSequenceRecord> ssrList = forheader.getFileHeader().getSequenceDictionary().getSequences();
		forheader.close();

		int chromcnt = 0;

		System.out.println("start realgin");
		int cnt = 0;
		for (SAMSequenceRecord ssr : ssrList) {

			String chrom = ssr.getSequenceName();
			boolean usualchrom = usualChrom(chrom);

			System.out.println("processing chr " + chrom);

			SAMFileReader normalbamr = getReader(normalbamf);
			SAMFileReader normalbamr2 = getReader(normalbamf);
			SAMFileReader tumorbamr = getReader(tumorbamf);
			SAMFileReader tumorbamr2 = getReader(tumorbamf);
			if (usualchrom) {

				// debug
				// if (!chrom.equals("chr17")) {
				// continue;
				// }

				realignChrom(chrom, normalbamr, normalbamw, tumorbamr, tumorbamw, tgr, dataset, numthread, normalbamr2,
						tumorbamr2);
				//
				// break;

			} else {

				copy(chrom, normalbamr, normalbamw);
				copy(chrom, tumorbamr, tumorbamw);

			}
			normalbamr.close();
			tumorbamr.close();
			normalbamr2.close();
			tumorbamr2.close();

		}

		normalbamw.close();
		tumorbamw.close();

	}

	public static final int FlgNormal = 1;
	public static final int FlgTumor = 2;

	private void realignChrom(String chrom, SAMFileReader normalbamr, SAMFileWriter normalbamw, SAMFileReader tumorbamr,
			SAMFileWriter tumorbamw, TwoBitGenomeReader tgr, DataSet dataset, int numthread, SAMFileReader normalbamr2,
			SAMFileReader tumorbamr2) throws Exception {

		//
		TreeMap<Integer, Indel> indelpos = new TreeMap<Integer, Indel>();
		int startR = 0;
		int endR = 0;

		int start = Integer.MAX_VALUE;
		int end = 0;

		// startR = 7675894;
		// endR = 7676762;

		// stats indel pos
		System.out.println("stat indel pos " + chrom);

		CloseableIterator<SAMRecord> iteN = normalbamr.query(chrom, startR, endR, false);
		Indel indel = null;
		while (iteN.hasNext()) {

			SAMRecord sam = iteN.next();
			if (!sam.getReadUnmappedFlag()) {
				if (sam.getAlignmentStart() < start) {
					start = sam.getAlignmentStart();
				}
				if (sam.getAlignmentEnd() > end) {
					end = sam.getAlignmentEnd();
				}
			}

			if (indel != null) {

				boolean intercect = intercect(sam, indel);
				if (intercect) {
					indel.incTotal();
				}

			}

			if (contatinIndel(sam)) {
				int ip = indelpos(sam);

				if (dataset != null) {
					if (null != dataset.getCh().getCapInterval(chrom, ip)) {

						if (indelpos.containsKey(ip)) {

							// indelpos.put(ip, indelpos.get(ip) + 1);
							indel = indelpos.get(ip);
							indel.inc(true);

						} else {

							indel = new Indel(sam);
							indel.inc(true);
							indelpos.put(ip, indel);
						}

					}
				} else {

					if (indelpos.containsKey(ip)) {

						// indelpos.put(ip, indelpos.get(ip) + 1);
						indel = indelpos.get(ip);
						indel.inc(true);

					} else {

						indel = new Indel(sam);
						indel.inc(true);
						indelpos.put(ip, indel);
					}

				}
			}

		}
		iteN.close();

		indel = null;
		CloseableIterator<SAMRecord> iteT = tumorbamr.query(chrom, startR, endR, false);
		while (iteT.hasNext()) {

			SAMRecord sam = iteT.next();
			if (!sam.getReadUnmappedFlag()) {
				if (sam.getAlignmentStart() < start) {
					start = sam.getAlignmentStart();
				}
				if (sam.getAlignmentEnd() > end) {
					end = sam.getAlignmentEnd();
				}
			}

			if (indel != null) {

				boolean intercect = intercect(sam, indel);
				if (intercect) {
					indel.incTotal();
				}

			}

			if (contatinIndel(sam)) {
				int ip = indelpos(sam);

				if (dataset != null) {
					if (null != dataset.getCh().getCapInterval(chrom, ip)) {

						if (indelpos.containsKey(ip)) {

							// indelpos.put(ip, indelpos.get(ip) + 1);
							indel = indelpos.get(ip);
							indel.inc(false);

						} else {

							indel = new Indel(sam);
							indel.inc(false);
							indelpos.put(ip, indel);
						}

					}
				} else {

					if (indelpos.containsKey(ip)) {

						// indelpos.put(ip, indelpos.get(ip) + 1);
						indel = indelpos.get(ip);
						indel.inc(false);

					} else {

						indel = new Indel(sam);
						indel.inc(false);
						indelpos.put(ip, indel);
					}

				}

			}

		}
		iteT.close();

		Set<Integer> single = new HashSet<Integer>();
		int cnt = 0;
		double frequencythres = 0.03;
		for (Entry<Integer, Indel> et : indelpos.entrySet()) {

			if (et.getValue().getCount() >= 5000) {
				single.add(et.getKey());
			}
			if (et.getValue().getCount() <= 2) {
				single.add(et.getKey());
			}
			if(et.getValue().getAlelleFrequency() <= frequencythres){
				
				//System.out.println("AF=" + et.getValue().getAlelleFrequency());
				single.add(et.getKey());
			}
			

		}
		for (int key : single) {
			indelpos.remove(key);
		}

		List<SAMRecord> normal = new ArrayList<SAMRecord>();
		List<SAMRecord> normal_realgin = new ArrayList<SAMRecord>();

		List<SAMRecord> tumor = new ArrayList<SAMRecord>();
		List<SAMRecord> tumor_realign = new ArrayList<SAMRecord>();
		//
		List<SAMRecord> realgin = new ArrayList<SAMRecord>();

		//
		// to Normal, Tumor, mixFor realgin
		System.out.println("extract realgin reads " + chrom);

		CloseableIterator<SAMRecord> iteN2 = normalbamr2.query(chrom, startR, endR, false);
		while (iteN2.hasNext()) {

			SAMRecord sam = iteN2.next();
			int nm = getNM(sam);
			// if depth is deep and indel is clear just take 1/3 of reads
			// to make it light weight
			if (contatinIndel(sam)) {

				if (enoughMatch(sam) && enoughDepth(sam, indelpos)) {

					//
					if (Math.random() < 0.85) {
						normal.add(sam);
						continue;
					}
				}

			}

			Indel idel = null;
			if (((nm > 0) && (contatinIndel(sam) || (idel = getIndel(sam, indelpos)) != null))) {

				sam.setAttribute("YY", FlgNormal);
				if (idel != null && idel.observedBoth()) {

					normal_realgin.add(sam);

				} else {

					realgin.add(sam);
				}

			} else {

				normal.add(sam);

			}

		}
		iteN2.close();

		CloseableIterator<SAMRecord> iteT2 = tumorbamr2.query(chrom, startR, endR, false);

		while (iteT2.hasNext()) {

			SAMRecord sam = iteT2.next();

			// if depth is deep and indel is clear just take 1/3 of reads
			// to make it light weight
			if (contatinIndel(sam)) {

				if (enoughMatch(sam) && enoughDepth(sam, indelpos)) {

					//
					if (Math.random() < 0.85) {
						tumor.add(sam);
						continue;
					}
				}

			}

			int nm = getNM(sam);

			Indel idel = null;
			if (((nm > 0) && (contatinIndel(sam) || (idel = getIndel(sam, indelpos)) != null))) {

				sam.setAttribute("YY", FlgTumor);
				if (idel != null && idel.observedBoth()) {

					tumor_realign.add(sam);

				} else {

					realgin.add(sam);

				}

			} else {

				tumor.add(sam);

			}

		}
		iteT2.close();

		if (realgin.size() > 1) {

			System.out.println("normal not to realgin " + normal.size());
			System.out.println("tumor not to realgin " + tumor.size());
			System.out.println("realign normal" + normal_realgin.size());
			System.out.println("realign tumor" + tumor_realign.size());
			System.out.println("realign normal with tumor" + realgin.size());

			System.out.println("sort reads for realgin " + chrom);
			// sort realgin
			Collections.sort(realgin, new SAMRecordCoordinateComparator());

			//

			System.out.println("realign by SRMA " + chrom + " size =" + realgin.size());

			TwobitReferenceSequence tbrs = new TwobitReferenceSequence(tgr, normalbamr.getFileHeader());

			ReferenceSequence res = tbrs.getSequence(chrom, start, end);

			// realgin dual re
			// ////////////////////////////////////
			realgin(chrom, normalbamr, numthread, normal, tumor, realgin, start, end, tbrs, res);

			// realgin dual re
			// ////////////////////////////////////
			realgin(chrom, normalbamr, numthread, normal, tumor, normal_realgin, start, end, tbrs, res);

			// realgin dual re
			// ////////////////////////////////////
			realgin(chrom, normalbamr, numthread, normal, tumor, tumor_realign, start, end, tbrs, res);

			System.out.println("sort normal " + chrom);
			// sort normal
			Collections.sort(normal, new SAMRecordCoordinateComparator());
			for (SAMRecord sam : normal) {

				normalbamw.addAlignment(sam);

			}
			normal = null;
			// normalbamw.close();

			System.out.println("sort tumor " + chrom);
			Collections.sort(tumor, new SAMRecordCoordinateComparator());
			for (SAMRecord sam : tumor) {

				tumorbamw.addAlignment(sam);

			}
			tumor = null;
			// tumorbamw.close();

		}

	}

	private boolean intercect(SAMRecord sam, Indel indel) {
		
		int start = sam.getAlignmentStart();
		int end = sam.getAlignmentEnd();
		
		int s = indel.pos;
		int e = indel.pos + indel.len;
		if(indel.insersion){
			e = indel.pos;
		}
		
		return s <=end && start<=e;
	}

	private void realgin(String chrom, SAMFileReader normalbamr, int numthread, List<SAMRecord> normal,
			List<SAMRecord> tumor, List<SAMRecord> realgin, int start, int end, TwobitReferenceSequence tbrs,
			ReferenceSequence res) throws InterruptedException {
		// sep list
		List<List<SAMRecord>> sep = sep(realgin);
		List<List<SAMRecord>> ret = new ArrayList<List<SAMRecord>>();
		ExecutorService exec = Executors.newFixedThreadPool(numthread);

		try {

			for (List<SAMRecord> list : sep) {

				Runtask task = new Runtask(list, ret, normalbamr.getFileHeader(), res, tbrs, chrom, start, end);

				exec.execute(task);

			}

		} finally {
			exec.shutdown();
			exec.awaitTermination(5, TimeUnit.HOURS);
		}

		System.out.println("realign by SRMA finish " + chrom);

		// ///////////////////////////////////
		//
		for (List<SAMRecord> rlist : ret) {
			// add to each
			for (SAMRecord sam : rlist) {

				int flg = sam.getIntegerAttribute("YY");
				if (flg == FlgNormal) {
					normal.add(sam);
				} else {
					tumor.add(sam);
				}

			}
		}
	}

	private boolean enoughDepth(SAMRecord sam, TreeMap<Integer, Indel> indelpos) {

		int s = sam.getAlignmentStart();
		int e = sam.getAlignmentEnd();
		//

		Integer f = indelpos.floorKey(e);
		Indel indel = null;
		if (f != null && (s < f && f < e)) {

			indel = indelpos.get(f);
		}

		Integer c = indelpos.ceilingKey(s);
		if (c != null && (s < c && c < e)) {

			indel = indelpos.get(f);
		}
		if (indel == null) {
			return false;
		}
		return indel.getCount() >= 100;
	}

	private boolean enoughMatch(SAMRecord sam) {

		List<CigarElement> l = sam.getCigar().getCigarElements();

		//
		if (l.size() < 3) {
			return false;
		}

		for (int n = 0; n < l.size() - 3; n++) {

			CigarElement c0 = l.get(n);
			CigarElement c1 = l.get(n + 1);
			CigarElement c2 = l.get(n + 2);

			if (c0.getOperator().equals(CigarOperator.M)
					&& (c1.getOperator().equals(CigarOperator.I) || c1.getOperator().equals(CigarOperator.D))
					&& c2.getOperator().equals(CigarOperator.M)) {

				if (c0.getLength() >= 30 || c2.getLength() >= 30) {

					return true;

				}

			}

		}
		return false;
	}

	private int getNM(SAMRecord sam) {
		Integer nm = sam.getIntegerAttribute("NM");

		if (nm != null) {
			return nm;
		}
		return 0;
	}

	class Runtask implements Runnable {

		List<SAMRecord> list;
		List<List<SAMRecord>> ret;
		SAMFileHeader fileHeader;
		ReferenceSequence rsf;
		ReferenceSequenceFile referenceSequenceFile;
		String chr;
		int start;
		int end;

		Runtask(List<SAMRecord> list, List<List<SAMRecord>> ret, SAMFileHeader fileHeader, ReferenceSequence rsf,
				ReferenceSequenceFile referenceSequenceFile, String chr, int start, int end) {

			this.list = list;
			this.ret = ret;
			this.fileHeader = fileHeader;
			this.rsf = rsf;
			this.referenceSequenceFile = referenceSequenceFile;
			this.chr = chr;
			this.start = start;
			this.end = end;

		}

		@Override
		public void run() {

			SRMA_mod inst = new SRMA_mod();
			try {
				List<SAMRecord> rlist = inst.doWorkSRMA(fileHeader, rsf, referenceSequenceFile, chr, start, end, list,
						1);

				synchronized (ret) {
					ret.add(rlist);
				}

			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

	public static int gap = 10000;

	private List<List<SAMRecord>> sep(List<SAMRecord> realgin) {

		List<List<SAMRecord>> l = new ArrayList<List<SAMRecord>>();

		int b4 = 0;
		List<SAMRecord> holder = new ArrayList<SAMRecord>();
		for (SAMRecord sam : realgin) {

			if (b4 == 0) {
				b4 = sam.getAlignmentEnd();
				continue;
			}

			if ((sam.getAlignmentStart() - b4) > gap) {

				//
				l.add(holder);
				holder = new ArrayList<SAMRecord>();
				holder.add(sam);

			} else {

				holder.add(sam);

			}

			b4 = sam.getAlignmentEnd();

		}
		l.add(holder);

		return l;
	}

	private Indel getIndel(SAMRecord sam, TreeMap<Integer, Indel> indelpos) {

		int s = sam.getAlignmentStart();
		int e = sam.getAlignmentEnd();
		//
		int readlen = sam.getReadLength();

		Integer c = indelpos.ceilingKey(s);
		Integer f = indelpos.floorKey(e);

		if (f != null) {
			if (Math.abs(s - f) < readlen) {

				Indel idel = indelpos.get(f);
				int sIn = idel.pos;
				int eIn = idel.pos + idel.len;
				if (idel.insersion) {
					sIn = idel.pos - idel.len;
				}

				//
				try {
					Set<Integer> mispos = getMisPos(sam);
					for (int n : mispos) {
						if (n >= sIn && n <= eIn) {
							return idel;
						}
					}
					//
					// return true;
				} catch (Exception ex) {
					ex.printStackTrace();
				}
			}
		}
		if (c != null) {
			if (Math.abs(c - e) < readlen) {

				Indel idel = indelpos.get(c);
				int sIn = idel.pos;
				int eIn = idel.pos + idel.len;
				try {
					Set<Integer> mispos = getMisPos(sam);
					for (int n : mispos) {
						if (n >= sIn && n <= eIn) {
							return idel;
						}
					}
					//
					return null;
				} catch (Exception ex) {
					ex.printStackTrace();
				}
				//
				// return true;
			}
		}
		return null;

	}

	private boolean nearIndel(SAMRecord sam, TreeMap<Integer, Indel> indelpos) {

		String readname2 = "HWI-D00677:85:CA2PBANXX:7:1115:18706:42919";
		if (sam.getReadName().equals(readname2)) {

			System.out.println("debug2 realgin");
			System.out.println("before=" + sam.getAlignmentStart() + " " + sam.getCigarString());

		}

		int s = sam.getAlignmentStart();
		int e = sam.getAlignmentEnd();
		//
		int readlen = sam.getReadLength();

		Integer c = indelpos.ceilingKey(s);
		Integer f = indelpos.floorKey(e);

		if (f != null) {
			if (Math.abs(s - f) < readlen) {

				Indel idel = indelpos.get(f);
				int sIn = idel.pos;
				int eIn = idel.pos + idel.len;
				if (idel.insersion) {
					sIn = idel.pos - idel.len;
				}

				//
				try {
					Set<Integer> mispos = getMisPos(sam);
					for (int n : mispos) {
						if (n >= sIn && n <= eIn) {
							return true;
						}
					}
					//
					// return true;
				} catch (Exception ex) {
					ex.printStackTrace();
				}
			}
		}
		if (c != null) {
			if (Math.abs(c - e) < readlen) {

				Indel idel = indelpos.get(c);
				int sIn = idel.pos;
				int eIn = idel.pos + idel.len;
				try {
					Set<Integer> mispos = getMisPos(sam);
					for (int n : mispos) {
						if (n >= sIn && n <= eIn) {
							return true;
						}
					}
					//
					return false;
				} catch (Exception ex) {
					ex.printStackTrace();
				}
				//
				// return true;
			}
		}
		return false;
	}

	/*
	 * Regexp for MD string.
	 *
	 * \G = end of previous match. (?:[0-9]+) non-capturing (why non-capturing?)
	 * group of digits. For this number of bases read matches reference. - or -
	 * Single reference base for case in which reference differs from read. - or
	 * - ^one or more reference bases that are deleted in read.
	 *
	 */
	static final Pattern mdPat = Pattern.compile("\\G(?:([0-9]+)|([ACTGNactgn])|(\\^[ACTGNactgn]+))");

	/**
	 * from picard code Produce reference bases from an aligned SAMRecord with
	 * MD string and Cigar.
	 * 
	 * @param rec
	 *            Must contain non-empty CIGAR and MD attribute.
	 * @param includeReferenceBasesForDeletions
	 *            If true, include reference bases that are deleted in the read.
	 *            This will make the returned array not line up with the read if
	 *            there are deletions.
	 * @return References bases corresponding to the read. If there is an
	 *         insertion in the read, reference contains '-'. If the read is
	 *         soft-clipped, reference contains '0'. If there is a skipped
	 *         region and includeReferenceBasesForDeletions==true, reference
	 *         will have Ns for the skipped region.
	 */
	public static Set<Integer> getMisPos(final SAMRecord rec) {

		Set<Integer> mispos = new HashSet<Integer>();
		int start = rec.getAlignmentStart();
		boolean includeReferenceBasesForDeletions = false;
		final String md = rec.getStringAttribute(SAMTag.MD.name());
		if (md == null) {
			throw new SAMException("Cannot create reference from SAMRecord with no MD tag, read: " + rec.getReadName());
		}
		// Not sure how long output will be, but it will be no longer than this.
		int maxOutputLength = 0;
		final Cigar cigar = rec.getCigar();
		if (cigar == null) {
			throw new SAMException("Cannot create reference from SAMRecord with no CIGAR, read: " + rec.getReadName());
		}
		for (final CigarElement cigarElement : cigar.getCigarElements()) {
			maxOutputLength += cigarElement.getLength();
		}
		final byte[] ret = new byte[maxOutputLength];
		int outIndex = 0;

		Matcher match = mdPat.matcher(md);
		int curSeqPos = 0;

		int savedBases = 0;
		final byte[] seq = rec.getReadBases();
		for (final CigarElement cigEl : cigar.getCigarElements()) {
			int cigElLen = cigEl.getLength();
			CigarOperator cigElOp = cigEl.getOperator();

			if (cigElOp == CigarOperator.SKIPPED_REGION) {

			}
			// If it consumes reference bases, it's either a match or a deletion
			// in the sequence
			// read. Either way, we're going to need to parse through the MD.
			else if (cigElOp.consumesReferenceBases()) {
				// We have a match region, go through the MD
				int basesMatched = 0;

				// Do we have any saved matched bases?
				while ((savedBases > 0) && (basesMatched < cigElLen)) {
					ret[outIndex++] = seq[curSeqPos++];
					savedBases--;
					basesMatched++;
				}

				while (basesMatched < cigElLen) {
					boolean matched = match.find();
					if (matched) {
						String mg;
						if (((mg = match.group(1)) != null) && (mg.length() > 0)) {
							// It's a number , meaning a series of matches
							int num = Integer.parseInt(mg);
							for (int i = 0; i < num; i++) {
								if (basesMatched < cigElLen) {
									ret[outIndex++] = seq[curSeqPos++];
								} else {
									savedBases++;
								}
								basesMatched++;
							}
						}

						else if (((mg = match.group(2)) != null) && (mg.length() > 0)) {
							// It's a single nucleotide, meaning a mismatch
							if (basesMatched < cigElLen) {
								ret[outIndex++] = StringUtil.charToByte(mg.charAt(0));
								curSeqPos++;
								mispos.add((start + outIndex));

							} else {
								// throw new IllegalStateException("Should never
								// happen.");
							}
							basesMatched++;
						} else if (((mg = match.group(3)) != null) && (mg.length() > 0)) {
							// It's a deletion, starting with a caret
							// don't include caret
							if (includeReferenceBasesForDeletions) {
								final byte[] deletedBases = StringUtil.stringToBytes(mg);
								System.arraycopy(deletedBases, 1, ret, outIndex, deletedBases.length - 1);
								outIndex += deletedBases.length - 1;
							}
							basesMatched += mg.length() - 1;

							// Check just to make sure.
							if (basesMatched != cigElLen) {
								// throw new SAMException("Got a deletion in
								// CIGAR (" + cigar + ", deletion " + cigElLen +
								// " length) with an unequal ref insertion in MD
								// (" + md + ", md " + basesMatched + "
								// length");
							}
							if (cigElOp != CigarOperator.DELETION) {
								// throw new SAMException ("Got an insertion in
								// MD ("+md+") without a corresponding deletion
								// in cigar ("+cigar+")");
							}

						} else {
							matched = false;
						}
					}

					if (!matched) {
						// throw new SAMException("Illegal MD pattern: " + md +
						// " for read " + rec.getReadName() +
						// " with CIGAR " + rec.getCigarString());
					}
				}

			} else if (cigElOp.consumesReadBases()) {
				// We have an insertion in read
				for (int i = 0; i < cigElLen; i++) {
					char c = (cigElOp == CigarOperator.SOFT_CLIP) ? '0' : '-';
					ret[outIndex++] = StringUtil.charToByte(c);
					curSeqPos++;
				}
			} else {
				// It's an op that consumes neither read nor reference bases. Do
				// we just ignore??
			}

		}
		return mispos;

		// if (outIndex < ret.length) {
		// byte[] shorter = new byte[outIndex];
		// System.arraycopy(ret, 0, shorter, 0, outIndex);
		// //return shorter;
		// }
		// return ret;
	}
	////

	private Integer indelpos(SAMRecord sam) {

		int pos = sam.getAlignmentStart();
		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.D)) {
				return pos;
			}
			if (ce.getOperator().equals(CigarOperator.I)) {
				return pos;
			}

			if (ce.getOperator().consumesReferenceBases()) {
				pos = pos + ce.getLength();
			}
		}

		return pos;

	}

	private boolean contatinIndel(SAMRecord sam) {

		//
		//

		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.D)) {
				return true;
			}
			if (ce.getOperator().equals(CigarOperator.I)) {
				return true;
			}

		}

		return false;

	}

	private void copy(String chrom, SAMFileReader bamr, SAMFileWriter bamw) {

		CloseableIterator<SAMRecord> ite = bamr.query(chrom, 0, 0, false);
		while (ite.hasNext()) {

			SAMRecord sam = ite.next();
			bamw.addAlignment(sam);

		}

	}

	private boolean usualChrom(String chrom) {

		String chromnum = chrom;
		if (chromnum.contains("chr")) {
			chromnum = chromnum.replaceAll("chr", "");
		}
		if (StringUtils.isNumeric(chromnum))
			return true;
		if (chromnum.equalsIgnoreCase("X"))
			return true;
		if (chromnum.equalsIgnoreCase("Y"))
			return true;

		return false;

	}

	public static Option getOption(String opt, String longOpt, boolean hasArg, String description, boolean required) {
		Option option = new Option(opt, longOpt, hasArg, description);
		option.setRequired(required);
		return option;
	}

	private static List<Option> getOptionListForKarkinos() {

		List<Option> optionlist = new ArrayList<Option>();
		optionlist.add(getOption("n", "normalBam", true, "normal bam file", true));
		optionlist.add(getOption("t", "tumorBam", true, "tumor bam file", true));

		optionlist.add(getOption("r", "reference", true, "2 bit genome reference file", true));

		optionlist.add(getOption("o", "outdir", true, "output directory", true));

		optionlist.add(getOption("ct", "captureTarget", true, "Capture target regions(bed format)", false));

		optionlist.add(getOption("nt", "num threads", true, "number of threads", false));

		return optionlist;

	}

}
