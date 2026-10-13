import java.util.HashMap;
import umich.ms.glyco.GlycanParser;
import umich.ms.glyco.GlycanResidue;

/**
 * FragPipe's MBG residue check before MBG runs (CmdMBGMatch.verifyResidues, CmdMBGMatch.java:213-229),
 * with the residue table GlycoMassLoader builds (GlycoMassLoader.java:64-77).
 *
 * Usage: java -cp "<fragpipe>/lib/*" MbgVerifyResidues.java glycan_residues.txt glycan_mods.txt "HexNAc(1),Hex(1)"
 * Exits 1 with FragPipe's message when a residue/mod string cannot be parsed as a glycan.
 */
public class MbgVerifyResidues {
    public static void main(String[] args) {
        HashMap<String, GlycanResidue> residues = GlycanParser.parseGlycoResiduesDB(args[0]);
        residues.putAll(GlycanParser.parseGlycoModsDB(args[1], residues.size(), residues));
        for (String s : args[2].split("[\\s,;/]+")) {
            if (GlycanParser.parseGlycanString(s, residues) == null) {
                System.err.printf("Error: MBG glycan difference input %s for expansion could not be parsed as a glycan. "
                    + "Please check that all residues/mods are defined and the format is correct.%n", s);
                System.exit(1);
            }
        }
    }
}
