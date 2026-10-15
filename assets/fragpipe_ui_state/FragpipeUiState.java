// FragpipeUiState: write FragPipe's saved UI state (fragpipe.workflow) for a run.
//
// FragPipe saves fragpipe.workflow before it runs any tool (FragpipeRun.java:371-375,
// saveRuntimeConfig -> MessageSaveUiState -> FragpipeCacheUtils.tabsSave): every UI key, with the
// UI default for keys the input workflow lacks, the run values (workdir, workflow.ram,
// workflow.threads, tool paths) and the tool-version header. This class drives FragPipe's own
// headless load path (Fragpipe.headless(), Fragpipe.java:473-534) and its save handler, so the
// file is FragPipe's, without running or validating any tool.
//
// Build (against the FragPipe 24.0 lib/ folder; the pipeline stages the .class):
//   javac --release 17 -cp "<fragpipe>/lib/*" FragpipeUiState.java
// Usage:
//   java -Djava.awt.headless=true -Duser.home=<dir> -cp "<fragpipe>/lib/*:<this dir>" FragpipeUiState \
//     <workflow> <manifest> <workdir> <tools-folder> <database or -> <out>
// workflow.ram / workflow.threads follow Fragpipe.headless() without --ram / --threads.
import java.nio.file.Paths;
import org.nesvilab.fragpipe.Fragpipe;
import org.nesvilab.fragpipe.FragpipeLoader;
import org.nesvilab.fragpipe.FragpipeLocations;
import org.nesvilab.fragpipe.api.Bus;
import org.nesvilab.fragpipe.api.PropsFile;
import org.nesvilab.fragpipe.messages.MessageLoadUi;
import org.nesvilab.fragpipe.messages.MessageManifestLoad;
import org.nesvilab.fragpipe.messages.MessageSaveUiState;
import org.nesvilab.fragpipe.tabs.TabWorkflow;

public class FragpipeUiState {
  public static void main(String[] a) throws Exception {
    if (a.length != 6) {
      System.err.println("Usage: FragpipeUiState <workflow> <manifest> <workdir> <tools-folder> <database or -> <out>");
      System.exit(2);
    }
    java.util.Locale.setDefault(java.util.Locale.US);  // FragPipeMain.main
    System.setProperty("java.awt.headless", "true");
    Fragpipe.headless = true;
    Fragpipe.workflowFile = Paths.get(a[0]);
    Fragpipe.manifestFile = Paths.get(a[1]);
    Fragpipe.workdir = Paths.get(a[2]).toAbsolutePath().normalize().toString();
    Fragpipe.toolsFolderPath = Paths.get(a[3]).toAbsolutePath().normalize().toString();

    // FragPipeMain.main resolves the Philosopher binary and version before startup
    java.lang.reflect.Method philosopher = Class.forName("org.nesvilab.fragpipe.FragPipeMain").getDeclaredMethod("getPhilosopherBin");
    philosopher.setAccessible(true);
    philosopher.invoke(null);

    // Fragpipe.main0 + Fragpipe.headless(), without MessageRun
    Bus.register(new FragpipeLoader());
    Fragpipe.initDone.await();
    PropsFile p = FragpipeLocations.get().tryLoadSilently(Fragpipe.workflowFile, "user");
    if (p.getProperty("workflow.ram") == null) {
      p.setProperty("workflow.ram", "0");
    }
    if (p.getProperty("workflow.threads") == null) {
      p.setProperty("workflow.threads", Math.max(1, Math.min(Runtime.getRuntime().availableProcessors() - 1, TabWorkflow.maxProcessors)) + "");
    }
    p.setProperty("workdir", Fragpipe.workdir);
    p.setProperty("fragpipe-config.tools-folder", Fragpipe.toolsFolderPath);
    if (!a[4].equals("-")) {
      p.setProperty("database.db-path", Paths.get(a[4]).toAbsolutePath().normalize().toString());
    }
    Bus.post(new MessageLoadUi(p, true, true));
    Bus.post(new MessageManifestLoad());
    Fragpipe.loadWorkflowDone.await();
    Fragpipe.loadManifestDone.await();
    Thread.sleep(500);  // as Fragpipe.headless()

    // FragpipeRun.saveRuntimeConfig: Fragpipe.on(MessageSaveUiState) runs on the posting thread
    Bus.post(new MessageSaveUiState(Paths.get(a[5])));
    System.exit(0);
  }
}
