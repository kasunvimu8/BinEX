import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.*;
import java.io.IOException;

public class PythonTestBed {
  private static final String PYTHON_PATH =
    "C:\\Users\\kasun\\AppData\\Local\\Programs\\Python\\Python38-32\\python ";
  private static final String SCRIPT_PATH = "\"C:\\Users\\kasun\\Desktop\\proposed-model\\mahalanobis_binning.py\" ";

  private Process executeProcess(List<String> arguments) {
    Process process = null;
    try {
      System.out.println(PYTHON_PATH + SCRIPT_PATH + getArguments(arguments));
      process = Runtime.getRuntime().exec(
        PYTHON_PATH + SCRIPT_PATH + getArguments(arguments),
        null);
      System.out.println(process);
    } catch (IOException e) {
      e.printStackTrace();
    }
    return process;
  }

  private String getArguments(List<String> arguments) {
    String space = " ";
    String quoteMark = "\"";
    String line = "";
    for (String argument : arguments) {
      line = line + space + quoteMark + argument + quoteMark;
    }
    return line;
  }

  public static void main(String[] args) throws IOException {
    PythonTestBed pythonTestBed = new PythonTestBed();

    List<String> arguments = new LinkedList<>();
	
    arguments.add("10s");
    arguments.add("sample_data/10s/10s_binned_contigs_features.csv");
    arguments.add("sample_data/10s/10s_unbinned_contigs_features.csv");
    arguments.add("sample_data/10s/10s_taxon.csv");
	arguments.add("sample_data/10s/10s_final_ouput.csv");
	
	/* System.out.println(pythonTestBed.getArguments(arguments)); */
	/* System.exit(0); */
	
    Process process = pythonTestBed.executeProcess(arguments);

    BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));

    BufferedReader standardError = new BufferedReader(new InputStreamReader(process.getErrorStream()));

    System.out.println("Here is the standard output of the command:\n");
    String s = null;

    while ((s = reader.readLine()) != null) {
      System.out.println(s);
    }
    System.out.println("Here is the standard error of the command (if any):\n");
    while ((s = standardError.readLine()) != null) {
      System.out.println(s);
    }
  }
}
