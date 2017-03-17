/*
 * Copyright Hiroki Ueda

 *  This program is free software; you can redistribute it and/or modify it under
 *	the terms of the GNU General Public License as published by the Free Software
 *	Foundation; either version 2 of the License, or (at your option) any later
 *	version.
	
 *	This program is distributed in the hope that it will be useful, but WITHOUT
 *	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *	FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 *	details.
	
 *	You should have received a copy of the GNU General Public License along with
 *	this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 *	Place, Suite 330, Boston, MA 02111-1307 USA
 */


package jp.ac.utokyo.rcast.realignimpl;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.CloseableIterator;


public class NearAdjust {

	static String rname = "HWI-HSQ700714:273:C5T04ACXX:7:2209:4376:45504";
	// GCGCTCCAGCTCGCGGTGGTGGTGGTGGGGTCTCCCTGCGTCTGTCCCTTAGACGTAGTCCTTGCGGTCGTAGCCTGTGCCCAGGCTGGCTCCCGGGCCG

	public static void main(String[] arg) {

	
	}

	

	public static SAMFileReader getReader(File INPUT) {

		// boolean validextension =
		// INPUT.getName().endsWith("sam")||INPUT.getName().endsWith("bam");
		// if(!validextension){
		//
		// }
		SAMFileReader reader = new SAMFileReader(INPUT);
		SAMFileHeader sfh = reader.getFileHeader();
		if (sfh.getAttribute("SO") == null || sfh.getAttribute("SO").equals("sorted")) {
			sfh.setSortOrder(SortOrder.coordinate);
		}
		reader.setValidationStringency(ValidationStringency.SILENT);

		return reader;
	}

	public static void fixNMForVolt(List<SAMRecord> retlist) {

		for (SAMRecord sam : retlist) {

			String pg = sam.getStringAttribute("PG");
			if ((pg != null) && (pg.equals("Volt") && (sam.getCigarLength() > 1))) {

				int nmm = 0;
				Integer nm = sam.getIntegerAttribute("NM");
				if (nm != null) {
					nmm = nm;
				}
				int add = getIndelLength(sam);
				nmm = nmm + add;
				sam.setAttribute("NM", nmm);

			}
		}

	}

	public static void fixMDandNM(List<SAMRecord> retlist, ReferenceSequence res) {

		for (SAMRecord sam : retlist) {

			String pg = sam.getStringAttribute("PG");
			if ((pg != null) && (pg.equals("Volt") && (sam.getCigarLength() > 1))) {

				int nmm = 0;
				Integer nm = sam.getIntegerAttribute("NM");
				if (nm != null) {
					nmm = nm;
				}
				int add = getIndelLength(sam);
				nmm = nmm + add;
				sam.setAttribute("NM", nmm);

			} else {
				fixNMandMD(sam, res);
			}
		}

	}

	public static int fixNMandMD(SAMRecord sam, ReferenceSequence res) {

		setMD(res, sam, true);
		int nm = getNM(sam.getStringAttribute("MD"));
		int nmWIndel = nm - getIndelLength(sam);
		sam.setAttribute("NM", nm);
		return nmWIndel;
	}

	public static final int margin50 = 50;
	

	private static boolean containI(SAMRecord sam) {
		for (CigarElement ce : sam.getCigar().getCigarElements()) {
			if (ce.getOperator().equals(CigarOperator.I)) {
				return true;
			}
		}
		return false;
	}

	private static boolean nearIndel(int nmWindel, int[] indels, SAMRecord sam) {

		if (nmWindel == 0)
			return false;
		if (getIndelLength(sam) == 0) {

			return indels[0] > 0 || indels[1] > 0;

		}
		return false;
	}

	private static int[] findIndelPos(List<SAMRecord> list, int idx, SAMRecord sam2) {

		int buf = 30;
		int start = idx - buf;
		if (start < 0)
			start = 0;
		int end = idx + buf;
		if (end >= list.size()) {
			end = list.size() - 1;
		}

		int[] ia = new int[] { 0, 0, 0 };
		int insertstart = 0;
		int deleationStart = 0;
		int deleationsize = 1;
		int delcnt = 0;
		int inscnt = 0;
		boolean find = true;
		for (int n = start; n <= end; n++) {

			SAMRecord sam = list.get(n);
			if(sam==sam2)continue;
			// System.out.println(sam.format());
			if (countIndel(sam) > 0) {

				int alstart = sam.getAlignmentStart();
				int genomeIdx = 0;
				for (CigarElement ce : sam.getCigar().getCigarElements()) {

					if (ce.getOperator().equals(CigarOperator.D)) {
						deleationStart = alstart + genomeIdx;
						deleationsize = ce.getLength();
						delcnt++;
						if (delcnt > 3) {
							find = true;
							break;
						}

					} else if (ce.getOperator().consumesReferenceBases()) {
						genomeIdx = genomeIdx + ce.getLength();
					}

				}

			}

			if (find == true) {
				ia[0] = 0;
				ia[1] = deleationStart;
				ia[2] = deleationsize;
			}
		}
		return ia;
	}

	private static int nmWithoutIndel(SAMRecord sam) {

		//
		int nm = sam.getIntegerAttribute("NM");
		int indeln = getIndelLength(sam);
		return nm - indeln;

	}

	private static void setMD(ReferenceSequence res, SAMRecord sam, boolean force) {
		// add md
		String md = sam.getStringAttribute("MD");
		if (md == null || force) {

			md = getMD(sam, res);
			sam.setAttribute("MD", md);
			int nm = getNM(md);
			sam.setAttribute("NM", nm);

		}

	}

	private static int getNM(String md) {

		//
		int n = 0;
		for (char ch : md.toCharArray()) {

			if ((ch == 'A') || (ch == 'T') || (ch == 'G') || (ch == 'C')) {
				n++;
			}

		}
		return n;
	}

	
	private static boolean nearIndelbutNoIndel(int nmExceptIndel, int indellength, int[] indels) {

		if (nmExceptIndel == 0)
			return false;
		if (indellength == 0) {

			if (indels[0] > 0 || indels[1] > 0) {
				return true;
			}

		}
		return false;
	}

	private static void softClip(SAMRecord sam) {

		int longestMidx = -1;
		int idx = 0;
		int maxsize = 0;
		if(sam.getCigar()==null)return;
		
		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.M)) {

				if (ce.getLength() > maxsize) {
					maxsize = ce.getLength();
					longestMidx = idx;

				}

			}

			idx++;
		}
		if (longestMidx == -1) {
			// no M
			// something wrong
			return;
		}

		Cigar cg = new Cigar();
		idx = 0;
		int readidx = 0;
		int refidx = 0;
		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (idx == longestMidx) {

				int als = sam.getAlignmentStart();
				sam.setAlignmentStart(als + refidx);
				cg.add(new CigarElement(readidx, CigarOperator.S));
				cg.add(ce);
				int left = sam.getReadLength() - (readidx + ce.getLength());
				if (left > 0) {
					cg.add(new CigarElement(left, CigarOperator.S));
				}
				break;
			}

			if (ce.getOperator().consumesReferenceBases()) {
				refidx = refidx + ce.getLength();
			}
			if (ce.getOperator().consumesReadBases()) {
				readidx = readidx + ce.getLength();
			}

			idx++;
		}

	}

	private static int countIndel(SAMRecord sam) {

		if(sam.getCigar()==null)return 0;
		int n = 0;
		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.D) || ce.getOperator().equals(CigarOperator.I)) {
				n++;
			}

		}
		return n;

	}

	private static void clearIfpossible(SAMRecord sam) {

		//
		if (!sam.getMateUnmappedFlag()) {

			if (sam.getMateReferenceIndex() == sam.getReferenceIndex()) {

				sam.setReadUmappedFlag(true);
				sam.setReferenceIndex(-1);
				sam.setAlignmentStart(0);

			}

		}

	}

	// private static void adNM(SAMRecord sam) {
	// // fix NM for Volt
	// String pg = sam.getStringAttribute("PG");
	// if ((pg != null) && (pg.equals("Volt") && (sam.getCigarLength() > 1))) {
	//
	// int nmm = 0;
	// Integer nm = sam.getIntegerAttribute("NM");
	// if (nm != null) {
	// nmm = nm;
	// }
	// int add = getIndelLength(sam);
	// nmm = nmm + add;
	// sam.setAttribute("NM", nmm);
	//
	// }
	//
	// }

	private static boolean softClipCheck(SAMRecord sam, ReferenceSequence res, int nmExceptIndel) {

		if (sam.getReadUnmappedFlag())
			return false;
		// if(true) return false;//debug

		boolean processed = false;

		if (sam.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.M)
				&& sam.getCigar().getCigarElement(0).getLength() >= 20) {

			int start = sam.getAlignmentStart();
			if (sam.getReadLength() <= 20)
				return false;
			int index = 0;
			int mismatch = 0;
			int lastmisidx = 0;
			int firstmatch = -1;
			for (int n = start; n < start + 15; n++) {

				if (index >= sam.getReadLength()) {
					break;
				}
				char c1 = sam.getReadString().charAt(index);
				char ref = (char) res.getBases()[n - 1];
				char c2 = sam.getReadString().charAt(index + 1);
				char ref2 = (char) res.getBases()[n];
				char c3 = sam.getReadString().charAt(index + 2);
				char ref3 = (char) res.getBases()[n + 1];
				boolean match0 = (c1 == ref);
				boolean match1 = (c2 == ref2);
				boolean match2 = (c3 == ref3);
				int mismatchcnt = 0;
				if (!match0)
					mismatchcnt++;
				if (!match1)
					mismatchcnt++;
				if (!match2)
					mismatchcnt++;
				boolean allsame = (c1 == c2) && (c2 == c3);

				if (match0 && match1 && match2) {
					if (firstmatch == -1) {
						firstmatch = index;
						break;
					}

				} else {
					if (c1 != ref) {
						mismatch++;
						lastmisidx = index;
					}

				}
				if (mismatchcnt >= 2) {

				} else {

					char c4 = sam.getReadString().charAt(index + 3);
					char ref4 = (char) res.getBases()[n + 2];
					int samecnt = 0;
					if (match0) {
						if (c1 == c4)
							samecnt++;
					}
					if (match1) {
						if (c2 == c4)
							samecnt++;
					}
					if (match2) {
						if (c3 == c4)
							samecnt++;
					}

					if (c4 == ref4 && samecnt < 2) {
						break;
					}
				}
				index++;
				processed = true;
			}

			if (firstmatch > 0) {

				int firstMlen = sam.getCigar().getCigarElements().get(0).getLength();
				if ((lastmisidx + 1) <= firstMlen) {
					int startn = sam.getAlignmentStart() + lastmisidx + 1;
					sam.setAlignmentStart(startn);
					Cigar cg = new Cigar();
					int n = 0;
					for (CigarElement ce : sam.getCigar().getCigarElements()) {

						if (n == 0) {
							//
							if ((lastmisidx + 1) > 0) {
								CigarElement ce0 = new CigarElement((lastmisidx + 1), CigarOperator.S);
								cg.add(ce0);
							}
							//
							if ((firstMlen - (lastmisidx + 1)) > 0) {
								CigarElement ce01 = new CigarElement((firstMlen - (lastmisidx + 1)), CigarOperator.M);
								cg.add(ce01);
							}

						} else {

							cg.add(ce);
						}
						n++;
					}
					sam.setCigar(cg);
					// sam.setCigarString((lastmisidx + 1) + "S" +
					// (sam.getReadLength() - (lastmisidx + 1)) + "M");
				}
			}

		}

		if (sam.getCigar().getCigarElement(sam.getCigarLength() - 1).getOperator().equals(CigarOperator.M)
				&& sam.getCigar().getCigarElement(sam.getCigarLength() - 1).getLength() >= 20) {

			int end = sam.getAlignmentEnd();
			int mismatch = 0;
			//
			int indexn = sam.getReadLength() - 1;
			int lastmisidxn = 0;
			int firstmatch = -1;
			for (int n = end; n >= end - 15; n--) {

				char c1 = sam.getReadString().charAt(indexn);
				// char ref =
				// res.getGenoemicNucFromRealPos(sam.getReferenceIndex(), n);
				char ref = (char) res.getBases()[n - 1];
				char c2 = sam.getReadString().charAt(indexn - 1);
				// char ref2 =
				// res.getGenoemicNucFromRealPos(sam.getReferenceIndex(), n-1);
				char ref2 = (char) res.getBases()[n - 2];

				char c3 = sam.getReadString().charAt(indexn - 2);
				char ref3 = (char) res.getBases()[n - 3];
				boolean match0 = (c1 == ref);
				boolean match1 = (c2 == ref2);
				boolean match2 = (c3 == ref3);
				int mismatchcnt = 0;
				if (!match0)
					mismatchcnt++;
				if (!match1)
					mismatchcnt++;
				if (!match2)
					mismatchcnt++;
				// boolean allsame = (c1 == c2) && (c2 ==c3);

				if (match0 && match1 && match2) {
					if (firstmatch == -1) {
						firstmatch = indexn;
						break;
					}

				} else {
					if (c1 != ref) {
						mismatch++;
						lastmisidxn = indexn;
					}

				}

				if (mismatchcnt >= 2) {

				} else {

					char c4 = sam.getReadString().charAt(indexn - 3);
					char ref4 = (char) res.getBases()[n - 4];
					int samecnt = 0;
					if (match0) {
						if (c1 == c4)
							samecnt++;
					}
					if (match1) {
						if (c2 == c4)
							samecnt++;
					}
					if (match2) {
						if (c3 == c4)
							samecnt++;
					}

					if (c4 == ref4 && samecnt < 2) {
						break;
					}
				}

				indexn--;

			}

			int lastMlen = sam.getCigar().getCigarElement(sam.getCigarLength() - 1).getLength();
			if ((firstmatch < sam.getReadLength() - 1)) {

				int sl = sam.getReadLength() - lastmisidxn;
				if (lastMlen >= sl && sl > 0 && lastmisidxn != 0) {

					// sam.setCigarString((sam.getReadLength() - sl) + "M" + sl
					// + "S");

					Cigar cg = new Cigar();
					int n = 0;
					for (CigarElement ce : sam.getCigar().getCigarElements()) {

						if (n == sam.getCigarLength() - 1) {
							//
							if ((lastMlen - sl) > 0) {
								CigarElement ce01 = new CigarElement((lastMlen - sl), CigarOperator.M);
								cg.add(ce01);
							}
							if (sl > 0) {
								CigarElement ce0 = new CigarElement(sl, CigarOperator.S);
								cg.add(ce0);
							}
						} else {

							cg.add(ce);
						}
						n++;
					}
					sam.setCigar(cg);
					processed = true;
				}
			}

		}
		return processed;
	}

	private static int getIndelLength(SAMRecord sam) {

		if(sam.getCigar()==null){
			return 0;
		}
		int n = 0;
		for (CigarElement ce : sam.getCigar().getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.D) || ce.getOperator().equals(CigarOperator.I)) {
				n = n + ce.getLength();
			}

		}
		return n;

	}

	private static int checkStart(int startn,String read, ReferenceSequence res,Cigar cg) {
		
		int max = 0;
		int maxmatch = 0;
		for(int n = -1;n<2;n++){
			
			int match = getMatch(startn+n,read,res, cg);
			if(match>maxmatch){
				maxmatch = match;
				max = n;
			}
		}
		
		return startn+max;
	}

	private static int getMatch(int start,String read, ReferenceSequence res,Cigar cg) {
		
		int localidx = 0;
		int s = start;
		int matchcount =0;
		for (CigarElement ce : cg.getCigarElements()) {

			if (ce.getOperator().equals(CigarOperator.M)) {

				int max = ce.getLength();
				for (int m = 0; m < max; m++) {

					if (localidx >= read.length())
						break;
					//
					char r = read.charAt(localidx);
					// char g =
					// res.getGenoemicNucFromRealPos(sam.getReferenceIndex(),
					// s);
					char g = (char) res.getBases()[s - 1];

					if (r == g) {
						matchcount++;
					} 
					localidx++;
					s++;
				}

			} else if (ce.getOperator().equals(CigarOperator.S)) {

				
				localidx = localidx + ce.getLength();

			} else if (ce.getOperator().consumesReadBases()) {
				// insersion
				
				localidx = localidx + ce.getLength();

			} else if (ce.getOperator().consumesReferenceBases()) {

				
				s = s + ce.getLength();

			}

		}
		return matchcount;
	}

	

	private static void adde(Cigar cg, int match, int ins, int del) {
		if (match != 0) {
			cg.add(new CigarElement(match, CigarOperator.MATCH_OR_MISMATCH));
		}
		if (ins != 0) {
			if (!cg.isEmpty()) {
				cg.add(new CigarElement(ins, CigarOperator.INSERTION));
			}
		}
		if (del != 0) {
			if (!cg.isEmpty()) {
				cg.add(new CigarElement(del, CigarOperator.DELETION));
			}
		}
	}

	public static String getMD(SAMRecord sam, ReferenceSequence res) {

		String read = sam.getReadString();
		Cigar cg = sam.getCigar();

		try {

			StringBuffer sb = new StringBuffer();

			int s = sam.getAlignmentStart();

			int localidx = 0;
			int matchcount = 0;

			for (CigarElement ce : cg.getCigarElements()) {

				if (ce.getOperator().equals(CigarOperator.M)) {

					int max = ce.getLength();
					for (int m = 0; m < max; m++) {

						if (localidx >= read.length())
							break;
						//
						char r = read.charAt(localidx);
						// char g =
						// res.getGenoemicNucFromRealPos(sam.getReferenceIndex(),
						// s);
						char g = (char) res.getBases()[s - 1];

						if (r != g) {

							if (matchcount > 0) {
								sb.append(matchcount + "" + g);
							} else {
								sb.append(g);
							}
							matchcount = 0;
						} else {
							matchcount++;
						}

						localidx++;
						s++;
					}

				} else if (ce.getOperator().equals(CigarOperator.S)) {

					if (matchcount > 0) {
						sb.append(matchcount);
						matchcount = 0;
					}
					localidx = localidx + ce.getLength();

				} else if (ce.getOperator().consumesReadBases()) {
					// insersion
					matchcount = matchcount + ce.getLength();
					localidx = localidx + ce.getLength();

				} else if (ce.getOperator().consumesReferenceBases()) {

					if (matchcount > 0) {
						sb.append(matchcount);
						matchcount = 0;
					}
					// String delref =
					// res.getMtbh().getGenomicSeq(sam.getReferenceIndex(), s +
					// 1, s + ce.getLength(),
					// true);
					String delref = getSubstring(res, s + 1, s + ce.getLength() + 1);

					// del
					sb.append("^" + delref);

					s = s + ce.getLength();

				}

			}
			if (matchcount > 0) {
				sb.append(matchcount);
				matchcount = 0;
			}
			// System.out.println(sb.toString());
			return sb.toString();
		} catch (Exception ex) {

		}
		return "";
	}

	private static String getSubstring(ReferenceSequence res, int i, int j) {

		try {
			StringBuffer sb = new StringBuffer();
			for (int n = i - 1; n < j - 1; n++) {
				sb.append((char) res.getBases()[n]);
			}
			return sb.toString();
		} catch (Exception ex) {

		}
		return "";

	}

}
