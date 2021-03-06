/*
 * LICENSE to be determined
 */
package srma;



import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequence;

import java.io.PrintStream;
import java.util.Iterator;
import java.util.List;

public class Alignment {
    public static final int GAP = '-';
    int length;
    byte read[];
    byte reference[];
    int positions[]; // one for each read base
    int positionsIndex[]; // one for each read base, index into read/reference

    public Alignment(SAMRecord record, ReferenceSequence sequence) throws Exception
    {
        Cigar cigar;
        List<CigarElement> cigarElements;
        Iterator<CigarElement> iter;
        CigarElement cigarElement;
        int cigarElementLength;
        CigarOperator cigarElementOperator;
        int readIndex, referenceIndex;
        int i, index;
        byte referenceBases[];
        byte readBases[];
        int alignmentStart;
        int positionsLength;

        alignmentStart = record.getAlignmentStart();
        cigar = record.getCigar();
        readBases = record.getReadBases();

        // Get alignment length
        this.length = positionsLength = 0;
        cigarElements = cigar.getCigarElements();
        iter = cigarElements.iterator();
        while(iter.hasNext()) {
            cigarElement = iter.next();
            switch(cigarElement.getOperator()) {
                case M:
                case I:
                case EQ: // will EQ exist ?
                case X: // will X exist ? 
                    positionsLength += cigarElement.getLength();
                    this.length += cigarElement.getLength();
                    break;
                case D:
                    this.length += cigarElement.getLength();
                    break;
                default:
                    break;
            }
        }

        this.read = new byte[this.length];
        this.reference = new byte[this.length];
        this.positions = new int[positionsLength];
        this.positionsIndex = new int[positionsLength];

        // Get reference bases
        referenceIndex = record.getReferenceIndex();
        if(referenceIndex < 0) {
            throw new Exception("Reference index out of range: " + referenceIndex);
        }
        referenceBases = sequence.getBases();

        // Copy over alignment
        iter = cigarElements.iterator();
        index = readIndex = referenceIndex = 0;
        while(iter.hasNext()) {
            cigarElement = iter.next();
            cigarElementLength = cigarElement.getLength();
            cigarElementOperator = cigarElement.getOperator();
            for(i=0;i<cigarElementLength;i++) {
                // TODO: make sure the bases are upper case
                switch(cigarElementOperator) {
                    case M:
                    case EQ: // will EQ exist ?
                    case X: // will X exist ? 
                        this.reference[index] = referenceBases[alignmentStart - 1 + referenceIndex];
                        this.read[index] = readBases[readIndex]; 
                        referenceIndex++;
                        readIndex++;
                        index++;
                        break;
                    case D:
                        this.reference[index] = referenceBases[alignmentStart - 1 + referenceIndex];
                        this.read[index] = Alignment.GAP;
                        referenceIndex++;
                        index++;
                        break;
                    case I:
                        this.reference[index] = Alignment.GAP; 
                        this.read[index] = readBases[readIndex]; 
                        readIndex++;
                        index++;
                        break;
                    case S:
                        // Ignore soft-clipping
                        readIndex++;
                        break;
                    case H:
                        // Ignore hard-clipping
                        break;
                    case P:
                        throw new AlignmentException("Illegal Cigar Operator: " + cigarElementOperator + " (Padded alignments not supported).");
                    case N:
                        throw new AlignmentException("Illegal Cigar Operator: " + cigarElementOperator + " (Spliced alignments not supported)."); 
                    default:
                        throw new AlignmentException("Illegal Cigar Operator: " + cigarElementOperator + " (Not supported)");
                }
            }
        }

        // left justify
        this.leftJustify();

        // Set positions etc;
        index = readIndex = 0;
        referenceIndex = alignmentStart - 1;
        while(index < this.length) {
            if(0 == index || Alignment.GAP != this.reference[index-1]) { // previous not an ins
                referenceIndex++;
            }
            if(Alignment.GAP != read[index]) { // M/EQ/X/I
                this.positions[readIndex] = referenceIndex;
                this.positionsIndex[readIndex] = index;
                readIndex++;
            }
            index++;
        }
    }

    private void leftJustify()
    {
        // Left-justify alignment
        int i;
        int prevDel, prevIns, startDel, endDel, startIns, endIns;

        i = prevDel = prevIns = 0;
        startDel = endDel = startIns = endIns = -1;

        while(i<this.length) {
            assert (0 == prevIns || 0 == prevDel);

            if(Alignment.GAP == this.read[i]) {
                if(0 == prevDel) {
                    startDel = i;
                }
                prevDel = 1;
                endDel = i;
                prevIns = 0;
                startIns = -1;
                endIns = -1;
                i++;
            }
            else if(Alignment.GAP == this.reference[i]) {
                if(0 == prevIns) {
                    startIns = i;
                }
                prevIns = 1;
                endIns = i;
                prevDel = 0;
                startDel = -1;
                endDel = -1;
                i++;
            }
            else {
                if(1 == prevDel) {
                    assert (0 < startDel);
                    assert (startDel <= endDel);
                    startDel--;
                    while(0 <= startDel && // Bases remaining to examine 
                            this.read[startDel] != Alignment.GAP && // Hit another deletion 
                            this.reference[startDel] != Alignment.GAP && // Hit an insertion 
                            this.reference[startDel] == this.reference[endDel]) { // src ref base matches dest ref base 
                        assert (Alignment.GAP != this.reference[startDel]);
                        assert (Alignment.GAP != this.reference[endDel]);
                        assert (Alignment.GAP != this.read[startDel]);
                        assert (Alignment.GAP == this.read[endDel]);
                        this.read[endDel] = this.read[startDel];
                        this.read[startDel] = Alignment.GAP;
                        startDel--;
                        endDel--;
                            }
                    endDel++; // We decremented when we exited the loop 
                    i = endDel;
                    assert (Alignment.GAP != this.read[i]);
                    assert (Alignment.GAP != this.reference[i]);
                }
                else if(1 == prevIns) {
                    assert (startIns <= endIns);
                    startIns--;
                    while(0 <= startIns && // Bases remaining to examine 
                            this.read[startIns] != Alignment.GAP && // Hit another deletion 
                            this.reference[startIns] != Alignment.GAP && // Hit an insertion 
                            this.read[startIns] == this.read[endIns]) { // src this.read base matches dest this.read base 
                        assert (Alignment.GAP != this.read[startIns]);
                        assert (Alignment.GAP != this.read[endIns]);
                        assert (Alignment.GAP != this.reference[startIns]);
                        assert (Alignment.GAP == this.reference[endIns]);
                        this.reference[endIns] = this.reference[startIns];
                        this.reference[startIns] = Alignment.GAP;
                        startIns--;
                        endIns--;
                            }
                    endIns++; // We decremented when we exited the loop 
                    i = endIns;
                    assert (Alignment.GAP != this.read[i]);
                    assert (Alignment.GAP != this.reference[i]);
                }
                else {
                    i++;
                }
                prevDel = 0;
                prevIns = 0;
                startDel = -1;
                endDel = -1;
                startIns = -1;
                endIns = -1;
            }
        }
    }

    public void print(PrintStream out)
    {
        int i, j;

        for(i=0;i<this.length;i++) {
            out.print((char)this.reference[i]);
        }
        out.println("");
        for(i=0;i<this.length;i++) {
            out.print((char)this.read[i]);
        }
        out.println("");
        for(i=j=0;i<this.length;i++) {
            if(Alignment.GAP != this.read[i]) {
                if(0 < j) {
                    out.print(",");
                }
                out.print(this.positions[j]);
                j++;
            }
        }
        out.println("");
        for(i=j=0;i<this.length;i++) {
            if(Alignment.GAP != this.read[i]) {
                if(0 < j) {
                    out.print(",");
                }
                out.print(this.positionsIndex[j]);
                j++;
            }
        }
        out.println("");
    }

    public class AlignmentException extends Exception {
        private static final long serialVersionUID = 1;

        public AlignmentException()
        {
            super();
        }

        public AlignmentException(String s)
        {
            super(s);
        }
    }
}
