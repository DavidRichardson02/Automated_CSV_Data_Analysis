//  DataSetModeling.c
//  Automated_CSV_Data_Analysis
//  DavidRichardson02


#include "DataSetModeling.h"
#include "CommonDefinitions.h"
#include "GeneralUtilities.h"
#include "StringUtilities.h"
#include "FileUtilities.h"
#include "DataExtraction.h"
#include "DataAnalysis.h"






void gaussian_erf_plot(const char* filePathName)
{
	double xmin = 0.0;
	double xmax = 5.0;
	double step = 0.01;
	
	int n; // Number of samples
	double* x = get_uniform_samples(&n, xmin, xmax, step);  // Sample the x-axis uniformly
	double* f = compute_gaussian(n, x); // Compute the Gaussian samples
	double* F = gaussian_riemann_sum_integration(n, x, f); // Compute the integral samples using Riemann approximation
	
	
	char *directoryPathName = find_file_directory_path(filePathName);
	char *gaussian_erf_data = combine_strings(directoryPathName, "gaussian_erf_data.txt");
	
	// Write the results to a file for plotting in MATLAB here
	FILE* data_file = fopen(gaussian_erf_data, "w");
	if (data_file == NULL)
	{
		perror("\n\nError: Could not open the file for writing in 'gaussian_erf_plot'.\n");
		exit(1);
	}
	for (int i = 0; i < n; i++)
	{
		fprintf(data_file, "%.17g %.17g %.17g\n", x[i], f[i], F[i]);
	}
	fclose(data_file);
	printf("\n\nData written to gaussian_erf_data.txt\n");
	
	
	char *matlabFile = combine_strings(directoryPathName, "MATLAB_Guassian_Erf_Data_Plot.m");
	// Using the contents of the file that was written to previously, 'data.txt', create a ".m" matlab script that when run from a matlab application will plot the graphs.
	FILE* matlab_file = fopen(matlabFile, "w"); ///*footnote 5 in 'Utilities.h'*
	if (matlab_file == NULL)
	{
		perror("\n\nError: Could not open the MATLAB script file for writing in 'gaussian_erf_plot'.\n");
		exit(1);
	}
	
	fprintf(matlab_file, "data = load('gaussian_erf_data.txt');\n");
	fprintf(matlab_file, "x = data(:,1);\n");
	fprintf(matlab_file, "f = data(:,2);\n");
	fprintf(matlab_file, "F = data(:,3);\n\n");
	fprintf(matlab_file, "f1 = figure('Name', 'Gaussian and Error Function Plot');\n");
	fprintf(matlab_file, "plot(x, f, 'r', 'LineWidth', 1.5);\n\n");
	fprintf(matlab_file, "hold on;\n");
	fprintf(matlab_file, "plot(x, F, 'b', 'LineWidth', 1.5);\n\n");
	fprintf(matlab_file, "xlabel('x');\n");
	fprintf(matlab_file, "ytit = '$\\frac{f}{F}$';\n");
	fprintf(matlab_file, "ylabel(ytit,'interpreter','latex');\n\n");
	fprintf(matlab_file, "t = 'Gaussian $f =  \\frac{2}{\\sqrt{\\pi}}  e^{-x^{2}} $ and Error Function: $F = erf(x)$';\n\n");
	fprintf(matlab_file, "title(t,'interpreter','latex');\n\n");
	fprintf(matlab_file, "grid on;\n");
	fclose(matlab_file);
	printf("\n\nMATLAB script written to Plot_Problem_One_MATLAB.m\n");
	
	
	free(directoryPathName);
	free(gaussian_erf_data);
	free(x);
	free(f);
	free(F);
}

void matlab_plot(const char* dataFilePathName, char **matlabScript, const char *matlabPlotDirectory)
{
	/*-----------   Write Data Set Fields into Seperate Files   -----------*/
	//const char *matlabPlotDirectory = create_directory(dataFilePathName);
	printf("\n\nmatlabPlotDirectory: %s", matlabPlotDirectory);
	const char *matlabFilePath = combine_strings(matlabPlotDirectory, "_Histogram_Data_MATLAB_Plot.m");
	// Write the results to a text file for plotting in MATLAB
	FILE* matlab_file = fopen(matlabFilePath, "w+");
	write_file_contents(matlabFilePath, matlabScript);
	fclose(matlab_file);
	printf("\n\nMATLAB script written to MATLAB_Histogram_Data_Plot.m\n");
}
char ** write_histogram_to_matlab_file(Histogram histogram, const char *filePathName, const char *fileNameMatlab, const char *matlabPlotDirectory)
{
	/*-----------   Write Data Set Fields into Seperate Files   -----------*/
	char *matlabPlotDirectoryFilePath = allocate_memory_char_ptr(strlen(filePathName) + strlen(fileNameMatlab) + 1);
	sprintf(matlabPlotDirectoryFilePath, "_MATLAB_Histogram_Data_Plot_%s.txt", fileNameMatlab);
	const char *matlabFilePath = combine_strings(matlabPlotDirectory, matlabPlotDirectoryFilePath);
	
	FILE* histogram_data_file = fopen(matlabFilePath, "a+");
	fprintf(histogram_data_file, "f1 = figure('Name', 'Particle Data Histogram Plots');\n");
	fprintf(histogram_data_file, "histogram = [");
	for (int i = 0; i < histogram.num_bins; i++)
	{
		fprintf(histogram_data_file, "%.17g ", histogram.min_value + i * histogram.bin_width);
	}
	fprintf(histogram_data_file, "];\n");
	
	
	fprintf(histogram_data_file, "counts = [");
	for (int i = 0; i < histogram.num_bins; i++)
	{
		fprintf(histogram_data_file, "%d ", histogram.bins[i]);
	}
	fprintf(histogram_data_file, "];\n");
	fprintf(histogram_data_file, "bar(histogram);\n");
	fprintf(histogram_data_file, "xlabel('%s');\n", fileNameMatlab);
	fprintf(histogram_data_file, "ylabel('Frequency');\n");
	fprintf(histogram_data_file, "title('%s %s');\n", fileNameMatlab, "Histogram");
	fclose(histogram_data_file);
	int lineCount = count_file_lines(matlabFilePath, MAX_NUM_FILE_LINES);
	char **matlabScript = read_file_contents(matlabFilePath, lineCount);
	return matlabScript;
}









/**
 * generate_individual_matlab_scripts
 *
 * Given a results directory and a list of field names, this function creates a
 * MATLAB script for each field (e.g., <fieldName>_plot.m) that:
 * 1. Loads the histogram data (<fieldName>_histogram.txt).
 * 2. Plots a histogram.
 * 3. Loads analysis data (<fieldName>_analysis.txt) to get mean and std.
 * 4. Overlays a normal distribution curve.
 */
void generate_individual_matlab_scripts(const char *resultsDir, char **fieldNames, int numFields)
{
	for (int i = 0; i < numFields; i++)
	{
		const char *fieldName = fieldNames[i];
		
		char matlabScriptPath[4096];
		snprintf(matlabScriptPath, sizeof(matlabScriptPath), "%s/%s_plot.m", resultsDir, fieldName);
		
		FILE *mf = fopen(matlabScriptPath, "a+");
		if (!mf)
		{
			fprintf(stderr, "Error creating MATLAB plot script: %s\n", matlabScriptPath);
			continue;
		}
		
		// MATLAB code generation
		fprintf(mf, "%% MATLAB script for plotting field: %s\n", fieldName);
		fprintf(mf, "hist_data = load('%s_histogram.txt');\n", fieldName);
		fprintf(mf, "bin_centers = (hist_data(:,1) + hist_data(:,2))/2;\n");
		fprintf(mf, "counts = hist_data(:,3);\n");
		fprintf(mf, "f%d = figure('Name', '%s Histogram');\n", i, fieldName);
		fprintf(mf, "bar(bin_centers, counts, 1.0, 'FaceColor',[0.7 0.7 0.7]);\n");
		fprintf(mf, "title('%s Histogram');\n", fieldName);
		fprintf(mf, "xlabel('%s Values');\n", fieldName);
		fprintf(mf, "ylabel('Frequency');\n");
		fprintf(mf, "grid on;\n");
		fprintf(mf, "hold on;\n");
		
		// Load analysis file (Mean, Std Dev)
		fprintf(mf, "%% Load analysis to overlay normal distribution\n");
		fprintf(mf, "A = fileread('%s_full_analysis.txt');\n", fieldName);
		fprintf(mf, "tokens_mean = regexp(A, 'Mean:\\s*(\\S+)', 'tokens');\n");
		fprintf(mf, "meanVal = 0; if ~isempty(tokens_mean), meanVal = str2double(tokens_mean{1}{1}); end\n");
		fprintf(mf, "tokens_std = regexp(A, 'Std Dev:\\s*(\\S+)', 'tokens');\n");
		fprintf(mf, "stdVal = 1; if ~isempty(tokens_std), stdVal = str2double(tokens_std{1}{1}); end\n");
		
		// Overlay normal distribution
		fprintf(mf, "x = linspace(min(bin_centers), max(bin_centers), 100);\n");
		fprintf(mf, "y = (1/(stdVal*sqrt(2*pi))) * exp(-0.5*((x-meanVal)/stdVal).^2);\n");
		fprintf(mf, "max_counts = max(counts);\n");
		fprintf(mf, "scale_factor = max_counts / max(y);\n");
		fprintf(mf, "plot(x, y*scale_factor, 'r-', 'LineWidth', 1.5);\n");
		fprintf(mf, "legend('Histogram', 'Normal Fit');\n");
		
		fclose(mf);
	}
}
void generate_comprehensive_matlab_script(const char *analysisDir, char **fieldNames, int numFields)
{
	char compScriptPath[4096];
	snprintf(compScriptPath, sizeof(compScriptPath), "%s/comprehensive_plots.m", analysisDir);
	
	FILE *comp = fopen(compScriptPath, "w+");
	if (!comp)
	{
		fprintf(stderr, "Error creating comprehensive_plots.m\n");
		return;
	}
	
	fprintf(comp, "%% Comprehensive MATLAB script to plot all analyzed fields separately\n");
	fprintf(comp, "close all; clear;\n\n");
	
	// fields cell array
	fprintf(comp, "fields = {");
	for (int i = 0; i < numFields; i++)
	{
		fprintf(comp, "'%s'", fieldNames[i]);
		if (i < numFields - 1) fprintf(comp, ", ");
	}
	fprintf(comp, "};\n");
	
	
	
	
	
	fprintf(comp, "numFields = length(fields);\n");
	fprintf(comp, "for i = 1:numFields\n");
	fprintf(comp, "    fieldName = fields{i};\n");
	fprintf(comp, "    hist_data = load([fieldName '_histogram.txt']);\n");
	fprintf(comp, "    bin_centers = (hist_data(:,1) + hist_data(:,2))/2;\n");
	fprintf(comp, "    counts = hist_data(:,3);\n");
	fprintf(comp, "    f = figure('Name', [fieldName ' Histogram']);\n");
	fprintf(comp, "    bar(bin_centers, counts, 1.0, 'FaceColor',[0.7 0.7 0.7]);\n");
	fprintf(comp, "    title([fieldName ' Histogram']);\n");
	fprintf(comp, "    xlabel('Value'); ylabel('Count');\n");
	fprintf(comp, "    grid on;\n");
	
	
	
	fprintf(comp, "    %% Load analysis (Mean, Std)\n");
	fprintf(comp, "    A = fileread([fieldName '_full_analysis.txt']);\n");
	fprintf(comp, "    tokens_mean = regexp(A, 'Mean:\\s*(\\S+)', 'tokens');\n");
	fprintf(comp, "    meanVal = 0; if ~isempty(tokens_mean), meanVal = str2double(tokens_mean{1}{1}); end\n");
	fprintf(comp, "    tokens_std = regexp(A, 'Std Dev:\\s*(\\S+)', 'tokens');\n");
	fprintf(comp, "    stdVal = 1; if ~isempty(tokens_std), stdVal = str2double(tokens_std{1}{1}); end\n");
	fprintf(comp, "    x = linspace(min(bin_centers), max(bin_centers), 100);\n");
	fprintf(comp, "    y = (1/(stdVal*sqrt(2*pi))) * exp(-0.5*((x-meanVal)/stdVal).^2);\n");
	fprintf(comp, "    max_counts = max(counts);\n");
	fprintf(comp, "    scale_factor = max_counts / max(y);\n");
	fprintf(comp, "	   hold on;\n");
	// hold off %any plots entered after “hold off” will be plotted on a new graph
	fprintf(comp, "    plot(x, y*scale_factor, 'r-', 'LineWidth', 1.5);\n");
	fprintf(comp, "    legend('Histogram', 'Normal Fit');\n");
	
	fprintf(comp, "    %% Save figure for this field\n");
	fprintf(comp, "    print([fieldName '_all'],'-dpng');\n");
	
	fprintf(comp, "end\n");
	
	fclose(comp);
}























/**
 * gaussian_erf_plot
 *
 * Existing example of generating/plotting Gaussian data.
 * In many workflows, you may keep it as-is or adapt it
 * to new naming conventions.
 */
void gaussian_error_function_plot(const char* filePathName)
{
	double xmin = 0.0;
	double xmax = 5.0;
	double step = 0.01;
	
	int n; // Number of samples
	double* x = get_uniform_samples(&n, xmin, xmax, step);
	double* f = compute_gaussian(n, x);
	double* F = gaussian_riemann_sum_integration(n, x, f);
	
	// Example: store output in a new "Modeling" subdirectory, if desired
	// e.g. /path/to/analysis_dir/Modeling/gaussian_erf_data.txt
	char *parentDir = find_file_directory_path(filePathName);
	char *modelingDir = combine_strings(parentDir, "Modeling/");
	create_directory(modelingDir, ""); // Make sure it exists
	
	// 1) Write the data
	char *gaussianErfDataPath = combine_strings(modelingDir, "gaussian_erf_data.txt");
	FILE* data_file = fopen(gaussianErfDataPath, "w");
	if (!data_file) {
		perror("Error: Could not open the file for writing in 'gaussian_erf_plot'.\n");
		exit(1);
	}
	for (int i = 0; i < n; i++) {
		fprintf(data_file, "%.17g %.17g %.17g\n", x[i], f[i], F[i]);
	}
	fclose(data_file);
	printf("\nData written to: %s\n", gaussianErfDataPath);
	
	// 2) Write a MATLAB script in the same Modeling directory
	char *matlabFilePath = combine_strings(modelingDir, "MATLAB_Gaussian_Erf_Data_Plot.m");
	FILE* mf = fopen(matlabFilePath, "w");
	if (!mf) {
		perror("Error: Could not open MATLAB script file for writing in 'gaussian_erf_plot'.\n");
		exit(1);
	}
	fprintf(mf, "data = load('%s');\n", gaussianErfDataPath);
	fprintf(mf, "x = data(:,1);\n");
	fprintf(mf, "f = data(:,2);\n");
	fprintf(mf, "F = data(:,3);\n\n");
	fprintf(mf, "fig = figure('Name', 'Gaussian and ERF');\n");
	fprintf(mf, "plot(x, f, 'r', 'LineWidth', 1.5);\n");
	fprintf(mf, "hold on;\n");
	fprintf(mf, "plot(x, F, 'b', 'LineWidth', 1.5);\n");
	fprintf(mf, "xlabel('x'); ylabel('Value');\n");
	fprintf(mf, "title('Gaussian f and ERF integral');\n");
	fprintf(mf, "grid on;\n");
	fprintf(mf, "legend('Gaussian','ERF');\n");
	fclose(mf);
	
	// Cleanup
	free(parentDir);
	free(modelingDir);
	free(gaussianErfDataPath);
	free(matlabFilePath);
	free(x);
	free(f);
	free(F);
}

/**
 * generate_fieldwise_modeling_scripts
 *
 * Formerly generate_individual_matlab_scripts(...).
 * This version:
 * 1) Looks for “<fieldName>_histogram.txt” and “<fieldName>_stats.txt” in `analysisDir`.
 * 2) Creates a single MATLAB script `<fieldName>_model.m` for each field.
 * 3) Overlays normal distribution, uses AD statistic, etc.
 */
void generate_fieldwise_modeling_scripts(const char *analysisDir, char **fieldNames, int numFields)
{
	// Optionally create a subdirectory for scripts
	// e.g., “analysisDir/MATLAB_Scripts”
	char scriptsDir[4096];
	snprintf(scriptsDir, sizeof(scriptsDir), "%s/MATLAB_Scripts/", analysisDir);
	create_directory(scriptsDir, "");
	
	for (int i = 0; i < numFields; i++)
	{
		const char *fieldName = fieldNames[i];
		
		// 1) Script name
		char scriptPath[4096];
		snprintf(scriptPath, sizeof(scriptPath), "%s/%s_model.m", scriptsDir, fieldName);
		
		FILE *mf = fopen(scriptPath, "w");
		if (!mf)
		{
			fprintf(stderr, "Error creating MATLAB script: %s\n", scriptPath);
			continue;
		}
		
		// 2) We assume you have:
		//    - <fieldName>_histogram.txt   => [binStart binEnd count]
		//    - <fieldName>_stats.txt       => (Mean, Std, Skewness, etc.)
		// Adjust filenames to match your code
		fprintf(mf, "%% MATLAB script for modeling field: %s\n", fieldName);
		fprintf(mf, "hist_data = load('%s/%s_histogram.txt');\n", analysisDir, fieldName);
		fprintf(mf, "bin_centers = (hist_data(:,1) + hist_data(:,2))/2;\n");
		fprintf(mf, "counts = hist_data(:,3);\n");
		fprintf(mf, "figure('Name', '%s Histogram');\n", fieldName);
		fprintf(mf, "bar(bin_centers, counts, 1.0, 'FaceColor',[0.7 0.7 0.7]);\n");
		fprintf(mf, "title('%s Histogram');\n", fieldName);
		fprintf(mf, "xlabel('%s');\n", fieldName);
		fprintf(mf, "ylabel('Frequency');\n");
		fprintf(mf, "grid on;\n");
		fprintf(mf, "hold on;\n\n");
		
		// 3) Load stats => Mean, Std, Skewness, AD
		fprintf(mf, "statsTxt = fileread('%s/%s_stats.txt');\n", analysisDir, fieldName);
		fprintf(mf, "exprMean = 'Mean:\\s*(\\S+)';\n");
		fprintf(mf, "exprStd  = 'Standard Deviation:\\s*(\\S+)';\n");
		fprintf(mf, "exprSkew = 'Skewness:\\s*(\\S+)';\n");
		fprintf(mf, "exprAD   = 'Anderson-Darling A\\^2:\\s*(\\S+)'; %% e.g.\n\n");
		fprintf(mf, "meanVal = 0; stdVal = 1; skewVal = 0; adVal = 0;\n");
		fprintf(mf, "tmp = regexp(statsTxt, exprMean, 'tokens'); if ~isempty(tmp), meanVal = str2double(tmp{1}{1}); end\n");
		fprintf(mf, "tmp = regexp(statsTxt, exprStd,  'tokens'); if ~isempty(tmp), stdVal  = str2double(tmp{1}{1}); end\n");
		fprintf(mf, "tmp = regexp(statsTxt, exprSkew, 'tokens'); if ~isempty(tmp), skewVal = str2double(tmp{1}{1}); end\n");
		fprintf(mf, "tmp = regexp(statsTxt, exprAD,   'tokens'); if ~isempty(tmp), adVal   = str2double(tmp{1}{1}); end\n\n");
		
		// 4) Overlay normal distribution
		fprintf(mf, "x = linspace(min(bin_centers), max(bin_centers), 200);\n");
		fprintf(mf, "y = (1/(stdVal*sqrt(2*pi))) * exp(-0.5*((x-meanVal)/stdVal).^2);\n");
		fprintf(mf, "max_counts = max(counts);\n");
		fprintf(mf, "scale_factor = max_counts / max(y);\n");
		fprintf(mf, "plot(x, y*scale_factor, 'r-', 'LineWidth', 1.5);\n");
		fprintf(mf, "legend('Histogram','Normal Fit');\n\n");
		
		// 5) Print out summary in MATLAB figure if you like
		fprintf(mf, "disp('======================');\n");
		fprintf(mf, "disp(['Field: %s']);\n", fieldName);
		fprintf(mf, "disp(['Mean = ' num2str(meanVal)]);\n");
		fprintf(mf, "disp(['Std  = ' num2str(stdVal)]);\n");
		fprintf(mf, "disp(['Skew = ' num2str(skewVal)]);\n");
		fprintf(mf, "disp(['Anderson-Darling = ' num2str(adVal)]);\n");
		fprintf(mf, "disp('======================');\n");
		
		fclose(mf);
	}
}

/**
 * generate_comprehensive_modeling_script
 *
 * Formerly generate_comprehensive_matlab_script(...).
 * Creates one big “comprehensive_plots.m” script that loops
 * over all fields. It’s basically a multi-figure approach.
 */
void generate_comprehensive_modeling_script(const char *analysisDir, char **fieldNames, int numFields)
{
	char compScriptPath[4096];
	snprintf(compScriptPath, sizeof(compScriptPath), "%s/comprehensive_model_plots.m", analysisDir);
	
	FILE *comp = fopen(compScriptPath, "w+");
	if (!comp)
	{
		fprintf(stderr, "Error creating comprehensive_model_plots.m\n");
		return;
	}
	
	fprintf(comp, "%% Comprehensive MATLAB script to plot all fields\n");
	fprintf(comp, "close all; clear;\n\n");
	
	// Put all field names into a cell array
	fprintf(comp, "fields = {");
	for (int i = 0; i < numFields; i++)
	{
		fprintf(comp, "'%s'", fieldNames[i]);
		if (i < numFields - 1) fprintf(comp, ", ");
	}
	fprintf(comp, "};\n");
	fprintf(comp, "numFields = length(fields);\n");
	fprintf(comp, "for i = 1:numFields\n");
	fprintf(comp, "    fieldName = fields{i};\n\n");
	
	// Load histogram data
	fprintf(comp, "    hist_data = load([fieldName '_histogram.txt']);\n");
	fprintf(comp, "    bin_centers = (hist_data(:,1) + hist_data(:,2))/2;\n");
	fprintf(comp, "    counts = hist_data(:,3);\n\n");
	fprintf(comp, "    figure('Name', [fieldName ' Histogram']);\n");
	fprintf(comp, "    bar(bin_centers, counts, 1.0, 'FaceColor',[0.7 0.7 0.7]);\n");
	fprintf(comp, "    title([fieldName ' Histogram']);\n");
	fprintf(comp, "    xlabel('Value'); ylabel('Count');\n");
	fprintf(comp, "    grid on; hold on;\n\n");
	
	// Load stats
	fprintf(comp, "    statsTxt = fileread([fieldName '_stats.txt']);\n");
	fprintf(comp, "    exprMean = 'Mean:\\s*(\\S+)';\n");
	fprintf(comp, "    exprStd  = 'Standard Deviation:\\s*(\\S+)';\n");
	fprintf(comp, "    exprSkew = 'Skewness:\\s*(\\S+)';\n");
	fprintf(comp, "    exprAD   = 'Anderson-Darling A\\^2:\\s*(\\S+)';\n");
	fprintf(comp, "    meanVal = 0; stdVal = 1; skewVal = 0; adVal = 0;\n");
	fprintf(comp, "    tmp = regexp(statsTxt, exprMean, 'tokens'); if ~isempty(tmp), meanVal = str2double(tmp{1}{1}); end\n");
	fprintf(comp, "    tmp = regexp(statsTxt, exprStd,  'tokens'); if ~isempty(tmp), stdVal  = str2double(tmp{1}{1}); end\n");
	fprintf(comp, "    tmp = regexp(statsTxt, exprSkew, 'tokens'); if ~isempty(tmp), skewVal = str2double(tmp{1}{1}); end\n");
	fprintf(comp, "    tmp = regexp(statsTxt, exprAD,   'tokens'); if ~isempty(tmp), adVal   = str2double(tmp{1}{1}); end\n\n");
	
	// Normal overlay
	fprintf(comp, "    x = linspace(min(bin_centers), max(bin_centers), 100);\n");
	fprintf(comp, "    y = (1/(stdVal*sqrt(2*pi))) * exp(-0.5*((x-meanVal)/stdVal).^2);\n");
	fprintf(comp, "    max_counts = max(counts);\n");
	fprintf(comp, "    scale_factor = max_counts / max(y);\n");
	fprintf(comp, "    plot(x, y*scale_factor, 'r-', 'LineWidth', 1.5);\n");
	fprintf(comp, "    legend('Histogram','Normal Fit');\n\n");
	fprintf(comp, "end\n");
	
	fclose(comp);
	printf("\nCreated a comprehensive model script: %s\n", compScriptPath);
}



































































/**
 * generate_modeling_matlab_scripts
 *
 * Creates MATLAB scripts for each plottable field in the 'modelingDir' directory.
 * Each script does some combination of:
 *   1. Load the histogram data from <fieldName>_histogram.txt.
 *   2. Load the analysis from <fieldName>_analysis.txt (Mean, Std, Skewness, etc.).
 *   3. Plot the histogram plus an overlaid normal distribution.
 *   4. Optionally, more advanced "modeling" lines or expansions (Gamma fits, Exponential, etc.).
 */
void generate_modeling_matlab_scripts(const char *modelingDir, char **fieldNames, int numFields)
{
	for (int i = 0; i < numFields; i++)
	{
		const char *fieldName = fieldNames[i];
		if (fieldName == NULL || strlen(fieldName) == 0)
			continue;
		
		// Construct the script path: e.g. <modelingDir>/<fieldName>_modeling.m
		char matlabScriptPath[4096];
		snprintf(matlabScriptPath, sizeof(matlabScriptPath), "%s/%s_modeling.m", modelingDir, fieldName);
		
		FILE *mf = fopen(matlabScriptPath, "w+");
		if (!mf)
		{
			fprintf(stderr, "Error creating MATLAB modeling script: %s\n", matlabScriptPath);
			continue;
		}
		
		fprintf(mf, "%% MATLAB script for modeling field: %s\n", fieldName);
		fprintf(mf, "%% This script loads histogram & analysis data, then plots histogram + normal dist.\n\n");
		
		// 1. Load histogram data  ( e.g. <fieldName>_histogram.txt => columns: [binStart, binEnd, count] )
		fprintf(mf, "hist_data = load(\'%s_histogram.txt\');\n", fieldName);
		fprintf(mf, "bin_starts = hist_data(:,1);\n");
		fprintf(mf, "bin_ends   = hist_data(:,2);\n");
		fprintf(mf, "counts     = hist_data(:,3);\n\n");
		fprintf(mf, "bin_centers = (bin_starts + bin_ends) / 2;\n");
		
		// 2. Load analysis data    ( e.g. <fieldName>_analysis.txt => lines with "Mean: X", "Std Dev: Y", etc. )
		fprintf(mf, "analysisText = fileread(\'%s_analysis.txt\');\n", fieldName);
		fprintf(mf, "exprMean = ’Mean:\\s*(\\\\S+)’;\n");
		fprintf(mf, "exprStd  = ’Std Dev:\\s*(\\\\S+)’;\n");
		fprintf(mf, "tokensMean = regexp(analysisText, exprMean, ’tokens’);\n");
		fprintf(mf, "tokensStd  = regexp(analysisText, exprStd,  ’tokens’);\n");
		fprintf(mf, "meanVal = 0; if ~isempty(tokensMean), meanVal = str2double(tokensMean{1}{1}); end\n");
		fprintf(mf, "stdVal  = 1; if ~isempty(tokensStd ), stdVal  = str2double(tokensStd{1}{1}); end\n\n");
		
		// 3. Plot the histogram
		fprintf(mf, "f = figure(\'Name\', \'%s Modeling\');\n", fieldName);
		fprintf(mf, "bar(bin_centers, counts, 1.0, \'FaceColor\', [0.7 0.7 0.7]);\n");
		fprintf(mf, "title(sprintf(\'%s Modeling: Histogram + Normal Fit\'));\n", fieldName);
		fprintf(mf, "xlabel(sprintf(\'%s Values\'));\n", fieldName);
		fprintf(mf, "ylabel(\'Count\');\n");
		fprintf(mf, "grid on;\n");
		fprintf(mf, "hold on;\n\n");
		
		// 4. Overlay a normal distribution
		//    we guess the data range from min(bin_starts) to max(bin_ends)
		fprintf(mf, "xMin = min(bin_starts);\n");
		fprintf(mf, "xMax = max(bin_ends);\n");
		fprintf(mf, "xVals = linspace(xMin, xMax, 200);\n");
		fprintf(mf, "yPdf  = (1.0/(stdVal * sqrt(2*pi))) * exp(-0.5*((xVals - meanVal)/stdVal).^2);\n");
		fprintf(mf, "%% scale the PDF to match histogram height\n");
		fprintf(mf, "maxCounts = max(counts);\n");
		fprintf(mf, "maxPdf    = max(yPdf);\n");
		fprintf(mf, "if maxPdf ~= 0\n");
		fprintf(mf, "   scaleFactor = maxCounts / maxPdf;\n");
		fprintf(mf, "else\n");
		fprintf(mf, "   scaleFactor = 1;\n");
		fprintf(mf, "end\n");
		fprintf(mf, "plot(xVals, yPdf * scaleFactor, \'r-\', \'LineWidth\', 1.5);\n");
		fprintf(mf, "legend(\'Histogram\', \'Normal Fit\', \'Location\', \'best\');\n");
		
		// 5. Optionally add extra “modeling” lines—for example, user-specified fits, param. lines, etc.
		// Example: a vertical line for the mean
		fprintf(mf, "plot([meanVal meanVal], [0 maxCounts], \'k--\', \'LineWidth\', 1.2);\n");
		fprintf(mf, "text(meanVal, maxCounts*0.9, sprintf(\'Mean=%.2f\', meanVal), \'HorizontalAlignment\', \'center\', \'Color\', \'k\');\n");
		
		// 6. Save .fig or .png automatically
		fprintf(mf, "%% Save the figure to <fieldName>_modeling.png\n");
		fprintf(mf, "saveas(f, sprintf(\'%s_modeling.png\'));\n", fieldName);
		
		fclose(mf);
	}
}






/**
 * generate_modeling_matlab_script_all_fields
 *
 * Creates a single all_fields_modeling.m script in the modelingDir
 * that plots each field’s histogram + normal overlay in a loop.
 */
void generate_modeling_matlab_script_all_fields(const char *modelingDir, char **fieldNames, int numFields)
{
	// e.g. .../Modeling/all_fields_modeling.m
	char compScriptPath[4096];
	snprintf(compScriptPath, sizeof(compScriptPath), "%s/all_fields_modeling.m", modelingDir);
	
	FILE *comp = fopen(compScriptPath, "w+");
	if (!comp)
	{
		fprintf(stderr, "Error creating all_fields_modeling.m\n");
		return;
	}
	
	fprintf(comp, "%% Comprehensive MATLAB script that loops over all plottable fields.\n");
	fprintf(comp, "close all; clearvars;\n\n");
	
	// Put the field names into a cell array in MATLAB:
	fprintf(comp, "fields = {\n");
	for (int i = 0; i < numFields; i++)
	{
		if (fieldNames[i] && strlen(fieldNames[i]) > 0) {
			fprintf(comp, "   \'%s\'\n", fieldNames[i]);
		}
	}
	fprintf(comp, "};\n\n");
	
	// Now the for-loop
	fprintf(comp, "for k = 1:length(fields)\n");
	fprintf(comp, "    fieldName = fields{k};\n");
	fprintf(comp, "    %% -- Load histogram\n");
	fprintf(comp, "    histData = load([fieldName \'_histogram.txt\']);\n");
	fprintf(comp, "    binStarts = histData(:,1);\n");
	fprintf(comp, "    binEnds   = histData(:,2);\n");
	fprintf(comp, "    counts    = histData(:,3);\n");
	fprintf(comp, "    binCenters = (binStarts + binEnds)/2;\n\n");
	
	fprintf(comp, "    %% -- Load analysis\n");
	fprintf(comp, "    txt = fileread([fieldName \'_analysis.txt\']);\n");
	fprintf(comp, "    reMean = \'Mean:\\\\s*(\\\\S+)\';\n");
	fprintf(comp, "    reStd  = \'Std Dev:\\\\s*(\\\\S+)\';\n");
	fprintf(comp, "    tokMean = regexp(txt, reMean, \'tokens\');\n");
	fprintf(comp, "    tokStd  = regexp(txt, reStd,  \'tokens\');\n");
	fprintf(comp, "    meanVal = 0; if ~isempty(tokMean), meanVal = str2double(tokMean{1}{1}); end\n");
	fprintf(comp, "    stdVal  = 1; if ~isempty(tokStd ), stdVal  = str2double(tokStd{1}{1}); end\n\n");
	
	fprintf(comp, "    figure(\'Name\', [fieldName \' Modeling\']);\n");
	fprintf(comp, "    bar(binCenters, counts, 1.0, \'FaceColor\', [0.6 0.6 0.9]);\n");
	fprintf(comp, "    title([fieldName \' - Hist + Normal Fit\']);\n");
	fprintf(comp, "    xlabel([fieldName \' Values\']); ylabel(\'Count\');\n");
	fprintf(comp, "    grid on; hold on;\n");
	fprintf(comp, "    xMin = min(binStarts); xMax = max(binEnds);\n");
	fprintf(comp, "    xVals = linspace(xMin, xMax, 300);\n");
	fprintf(comp, "    yPdf = (1./(stdVal.*sqrt(2*pi))) .* exp(-0.5.*((xVals - meanVal)./stdVal).^2);\n");
	fprintf(comp, "    maxC = max(counts); maxY = max(yPdf);\n");
	fprintf(comp, "    sc   = maxC/maxY;\n");
	fprintf(comp, "    plot(xVals, yPdf*sc, \'r-\', \'LineWidth\', 1.5);\n");
	fprintf(comp, "    legend(\'Histogram\', \'Normal Fit\');\n");
	fprintf(comp, "    %% Save figure:\n");
	fprintf(comp, "    saveas(gcf, [fieldName \'_all_fields_plot.png\']);\n\n");
	fprintf(comp, "end\n");
	
	fclose(comp);
}
