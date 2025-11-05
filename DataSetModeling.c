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








/**
 * emit_adv_modeling_snippet_individual
 *
 * Writes an optional advanced-modeling MATLAB segment into a per-field script. The
 * emitted code attempts to load raw samples, fit several candidate distributions
 * (Normal, Lognormal, Gamma, Exponential when valid), select the best by AIC, and
 * overlay the best-fit PDF on the histogram. It also generates ECDF and QQ plots and
 * saves them alongside the script. All operations are guarded to remain no-risk when
 * required toolboxes are missing or files are absent.
 *
 * Works by printing MATLAB commands via fprintf into the open script stream `mf`.
 * The MATLAB snippet resolves raw data from multiple candidate paths, computes
 * log-likelihoods per fit, ranks via AIC, overlays the scaled best PDF, updates the
 * legend, and conditionally saves ECDF/QQ figure files. Errors are caught and surfaced
 * as MATLAB warnings at runtime.
 *
 * @note The emitted MATLAB expects variables `analysisDir`, `scriptDir`, `field`,
 *       `bin_centers`, and `counts` to be defined earlier in the script.
 *
 * @param mf Pointer to an open FILE stream for the per-field MATLAB script.
 * @return void.
 */
static void emit_adv_modeling_snippet_individual(FILE *mf)
{
	fprintf(mf,
			"%% --- Optional advanced modeling (raw samples + toolbox, no-risk) ---\n"
			"enableAdvancedFits = true;\n"
			"rawCandidates = {\n"
			"  fullfile(analysisDir, sprintf('%%s_samples.txt', field)),\n"
			"  sprintf('%%s_samples.txt', field)\n"
			"};\n"
			"rawLoaded = false;\n"
			"for rk=1:numel(rawCandidates)\n"
			"  if exist(rawCandidates{rk},'file'), try raw = load(rawCandidates{rk}); rawLoaded=true; break; end, end\n"
			"end\n"
			"if enableAdvancedFits && rawLoaded && exist('fitdist','file') && exist('aicbic','file')\n"
			"  try\n"
			"    Cnames = {}; Cdist = {}; LL = []; K = [];  %% names, dists, loglik, k-params\n"
			"    %% Normal (k=2)\n"
			"    d = fitdist(raw(:),'Normal');      Cnames{end+1}='Normal';    Cdist{end+1}=d; LL(end+1)=sum(log(pdf(d,raw(:))+eps)); K(end+1)=2;\n"
			"    if all(raw(:)>0)\n"
			"      %% Lognormal (k=2)\n"
			"      d = fitdist(raw(:),'Lognormal'); Cnames{end+1}='Lognormal'; Cdist{end+1}=d; LL(end+1)=sum(log(pdf(d,raw(:))+eps)); K(end+1)=2;\n"
			"      %% Gamma (k=2)\n"
			"      d = fitdist(raw(:),'Gamma');     Cnames{end+1}='Gamma';     Cdist{end+1}=d; LL(end+1)=sum(log(pdf(d,raw(:))+eps)); K(end+1)=2;\n"
			"      %% Exponential (k=1)\n"
			"      d = fitdist(raw(:),'Exponential'); Cnames{end+1}='Exponential'; Cdist{end+1}=d; LL(end+1)=sum(log(pdf(d,raw(:))+eps)); K(end+1)=1;\n"
			"    end\n"
			"    [AIC,~] = aicbic(LL, K, numel(raw)); [~,ix] = min(AIC);\n"
			"    bestName = Cnames{ix}; bestDist = Cdist{ix};\n"
			"    xfine = linspace(min(bin_centers), max(bin_centers), 300);\n"
			"    ypdf  = pdf(bestDist, xfine); sf2 = max(counts) / max(ypdf + eps);\n"
			"    plot(xfine, ypdf*sf2, 'm-', 'LineWidth', 1.2);\n"
			"    legend('Histogram','Normal Fit','Mean',['Best: ' bestName],'Location','best');\n"
			"    %% ECDF + QQ (saved next to script; guards keep Octave OK)\n"
			"    f_ecdf = figure('Name',[field ' ECDF']); ecdf(raw(:)); grid on;\n"
			"    saveas(f_ecdf, fullfile(scriptDir, sprintf('%%s_ecdf.png', field)));\n"
			"    if exist('savefig','file'), savefig(f_ecdf, fullfile(scriptDir, sprintf('%%s_ecdf.fig', field))); end\n"
			"    f_qq = figure('Name',[field ' QQ Normal']); qqplot(raw(:)); grid on;\n"
			"    saveas(f_qq,  fullfile(scriptDir, sprintf('%%s_qq.png', field)));\n"
			"    if exist('savefig','file'), savefig(f_qq,  fullfile(scriptDir, sprintf('%%s_qq.fig', field))); end\n"
			"  catch ME\n"
			"    warning('Advanced modeling failed for %%s: %%s', field, ME.message);\n"
			"  end\n"
			"end\n");
}


/**
 * emit_adv_modeling_snippet_comprehensive
 *
 * Writes an optional advanced-modeling MATLAB segment into the comprehensive script
 * (multi-field loop). The emitted code mirrors the per-field version: it locates raw
 * samples, fits candidate distributions, selects the AIC winner, overlays its PDF,
 * and saves ECDF/QQ plots, all guarded for toolbox/file availability.
 *
 * Works by printing MATLAB code into `comp`. At runtime, for each field `fn`, the
 * snippet loads raw samples if present, evaluates candidate fits, computes AIC,
 * overlays the best PDF scaled to the histogram, and saves ECDF/QQ figures with
 * names derived from `fn`. Errors are caught and reported as MATLAB warnings.
 *
 * @note The emitted MATLAB expects `analysisDir`, `scriptDir`, `fn`, `bc`, and `cnt`
 *       to exist in the surrounding script loop.
 *
 * @param comp Pointer to an open FILE stream for the comprehensive MATLAB script.
 * @return void.
 */
static void emit_adv_modeling_snippet_comprehensive(FILE *comp)
{
	fprintf(comp,
			"%% --- Optional advanced modeling (raw samples + toolbox, no-risk) ---\n"
			"enableAdvancedFits = true;\n"
			"rawCandidates = {\n"
			"  fullfile(analysisDir, sprintf('%%s_samples.txt', fn)),\n"
			"  sprintf('%%s_samples.txt', fn)\n"
			"};\n"
			"rawLoaded = false;\n"
			"for rk=1:numel(rawCandidates)\n"
			"  if exist(rawCandidates{rk},'file'), try raw = load(rawCandidates{rk}); rawLoaded=true; break; end, end\n"
			"end\n"
			"if enableAdvancedFits && rawLoaded && exist('fitdist','file') && exist('aicbic','file')\n"
			"  try\n"
			"    Cnames = {}; Cdist = {}; LL = []; K = [];\n"
			"    d = fitdist(raw(:),'Normal');      Cnames{end+1}='Normal';    Cdist{end+1}=d; LL(end+1)=sum(log(pdf(d,raw(:))+eps)); K(end+1)=2;\n"
			"    if all(raw(:)>0)\n"
			"      d = fitdist(raw(:),'Lognormal'); Cnames{end+1}='Lognormal'; Cdist{end+1}=d; LL(end+1)=sum(log(pdf(d,raw(:))+eps)); K(end+1)=2;\n"
			"      d = fitdist(raw(:),'Gamma');     Cnames{end+1}='Gamma';     Cdist{end+1}=d; LL(end+1)=sum(log(pdf(d,raw(:))+eps)); K(end+1)=2;\n"
			"      d = fitdist(raw(:),'Exponential'); Cnames{end+1}='Exponential'; Cdist{end+1}=d; LL(end+1)=sum(log(pdf(d,raw(:))+eps)); K(end+1)=1;\n"
			"    end\n"
			"    [AIC,~] = aicbic(LL, K, numel(raw)); [~,ix] = min(AIC);\n"
			"    bestName=Cnames{ix}; bestDist=Cdist{ix};\n"
			"    xfine = linspace(min(bc), max(bc), 300);\n"
			"    ypdf  = pdf(bestDist, xfine); sf2 = max(cnt)/max(ypdf+eps);\n"
			"    plot(xfine, ypdf*sf2, 'm-', 'LineWidth', 1.2);\n"
			"    legend('Histogram','Normal Fit','Mean',['Best: ' bestName],'Location','best');\n"
			"    f_ecdf = figure('Name',[fn ' ECDF']); ecdf(raw(:)); grid on;\n"
			"    saveas(f_ecdf, fullfile(scriptDir, sprintf('%%s_ecdf.png', fn)));\n"
			"    if exist('savefig','file'), savefig(f_ecdf, fullfile(scriptDir, sprintf('%%s_ecdf.fig', fn))); end\n"
			"    f_qq = figure('Name',[fn ' QQ Normal']); qqplot(raw(:)); grid on;\n"
			"    saveas(f_qq,  fullfile(scriptDir, sprintf('%%s_qq.png', fn)));\n"
			"    if exist('savefig','file'), savefig(f_qq,  fullfile(scriptDir, sprintf('%%s_qq.fig', fn))); end\n"
			"  catch ME\n"
			"    warning('Advanced modeling failed for %%s: %%s', fn, ME.message);\n"
			"  end\n"
			"end\n");
}




/**
 * emit_kde_overlay_individual
 *
 * Writes an optional KDE overlay section to a per-field MATLAB script. If
 * `ksdensity` is available, it estimates a smooth density either from loaded raw
 * samples or from pseudo-samples expanded from the histogram bins, scales it to the
 * histogram height, and overlays the curve with an updated legend.
 *
 * Works by printing MATLAB commands that: attempt to load raw samples from multiple
 * locations; fall back to expanding bin centers by counts; compute KDE on a fine grid;
 * scale and plot the KDE; and guard the entire section with try/catch and
 * `exist('ksdensity','file')`.
 *
 * @note Requires `field`, `analysisDir`, `bin_centers`, and `counts` to be defined
 *       earlier in the script that includes this snippet.
 *
 * @param mf Pointer to an open FILE stream for the per-field MATLAB script.
 * @return void.
 */
static void emit_kde_overlay_individual(FILE *mf)
{
	fprintf(mf,
			"%%%% %% Section: KDE Overlay (optional)\n"
			"if exist('ksdensity','file')\n"
			"  try\n"
			"    %%%% If raw samples exist, prefer them; else approximate from histogram\n"
			"    haveRaw=false; rawCandidates={fullfile(analysisDir,sprintf('%%s_samples.txt',field)),sprintf('%%s_samples.txt',field)};\n"
			"    for rk=1:numel(rawCandidates)\n"
			"      if exist(rawCandidates{rk},'file'), try raw=load(rawCandidates{rk}); haveRaw=true; break; end, end\n"
			"    end\n"
			"    if haveRaw\n"
			"      xk = linspace(min(bin_centers), max(bin_centers), 400);\n"
			"      [yk,~] = ksdensity(raw(:), xk);\n"
			"    else\n"
			"      %%%% expand counts into pseudo-samples (coarse but OK if counts small)\n"
			"      ps=[]; for ii=1:numel(bin_centers), ps=[ps; repmat(bin_centers(ii), counts(ii), 1)]; end\n"
			"      if ~isempty(ps)\n"
			"        xk = linspace(min(bin_centers), max(bin_centers), 400);\n"
			"        [yk,~] = ksdensity(ps, xk);\n"
			"      else, xk=[]; yk=[]; end\n"
			"    end\n"
			"    if ~isempty(yk)\n"
			"      sfK = max(counts) / max(yk + eps);\n"
			"      plot(xk, yk*sfK, 'g-', 'LineWidth', 1.2);\n"
			"      legend('Histogram','Normal Fit','Mean','KDE','Location','best');\n"
			"    end\n"
			"  catch ME\n"
			"    warning('KDE overlay skipped: %s', ME.message);\n"
			"  end\n"
			"end\n\n"
			);
}


/**
 * emit_kde_overlay_comprehensive
 *
 * Writes an optional KDE overlay section to the comprehensive MATLAB script. Behavior
 * mirrors the per-field version but operates on variables for the current field in
 * the loop. If `ksdensity` is present, it overlays a scaled KDE curve and updates
 * the legend; otherwise the code path is skipped safely.
 *
 * Works by emitting guarded MATLAB code that prefers raw samples and falls back to
 * pseudo-samples expanded from `(bc, cnt)`, then plots a scaled KDE on top of the
 * histogram.
 *
 * @note Expects `fn`, `analysisDir`, `bc`, and `cnt` to be defined by the caller.
 *
 * @param comp Pointer to an open FILE stream for the comprehensive MATLAB script.
 * @return void.
 */
static void emit_kde_overlay_comprehensive(FILE *comp)
{
	fprintf(comp,
			"%%%% %% Section: KDE Overlay (optional)\n"
			"if exist('ksdensity','file')\n"
			"  try\n"
			"    haveRaw=false; rawCandidates={fullfile(analysisDir,sprintf('%%s_samples.txt',fn)),sprintf('%%s_samples.txt',fn)};\n"
			"    for rk=1:numel(rawCandidates)\n"
			"      if exist(rawCandidates{rk},'file'), try raw=load(rawCandidates{rk}); haveRaw=true; break; end, end\n"
			"    end\n"
			"    if haveRaw\n"
			"      xk = linspace(min(bc), max(bc), 400);\n"
			"      [yk,~] = ksdensity(raw(:), xk);\n"
			"    else\n"
			"      ps=[]; for ii=1:numel(bc), ps=[ps; repmat(bc(ii), cnt(ii), 1)]; end\n"
			"      if ~isempty(ps)\n"
			"        xk = linspace(min(bc), max(bc), 400);\n"
			"        [yk,~] = ksdensity(ps, xk);\n"
			"      else, xk=[]; yk=[]; end\n"
			"    end\n"
			"    if ~isempty(yk)\n"
			"      sfK = max(cnt) / max(yk + eps);\n"
			"      plot(xk, yk*sfK, 'g-', 'LineWidth', 1.2);\n"
			"      legend('Histogram','Normal Fit','Mean','KDE','Location','best');\n"
			"    end\n"
			"  catch ME\n"
			"    warning('KDE overlay skipped: %s', ME.message);\n"
			"  end\n"
			"end\n\n"
			);
}







/**
 * emit_gmm_overlay_individual
 *
 * Writes an optional Gaussian Mixture Model (1–3 components) overlay section for a
 * per-field MATLAB script. When `fitgmdist` is available and enough samples exist, it
 * fits mixtures with k ∈ {1,2,3}, selects the model with the lowest BIC, and overlays
 * its PDF scaled to the histogram, extending the legend accordingly.
 *
 * Works by printing MATLAB code that loads or synthesizes samples, iterates k, fits
 * mixtures with regularization, tracks BIC, and overlays the best model’s PDF. Errors
 * are caught and only a warning is emitted at runtime.
 *
 * @note Expects `analysisDir`, `field`, `bin_centers`, and `counts` to exist in the
 *       host MATLAB script.
 *
 * @param mf Pointer to an open FILE stream for the per-field MATLAB script.
 * @return void.
 */
static void emit_gmm_overlay_individual(FILE *mf)
{
	fprintf(mf,
			"%%%% %% Section: Gaussian Mixture Model (BIC, 1..3 comps)\n"
			"if exist('fitgmdist','file')\n"
			"  try\n"
			"    haveRaw=false; rawCandidates={fullfile(analysisDir,sprintf('%%s_samples.txt',field)),sprintf('%%s_samples.txt',field)};\n"
			"    for rk=1:numel(rawCandidates)\n"
			"      if exist(rawCandidates{rk},'file'), try raw=load(rawCandidates{rk}); haveRaw=true; break; end, end\n"
			"    end\n"
			"    if ~haveRaw\n"
			"      raw=[]; for ii=1:numel(bin_centers), raw=[raw; repmat(bin_centers(ii), counts(ii), 1)]; end\n"
			"    end\n"
			"    if numel(raw)>10\n"
			"      bestGM=[]; bestBIC=Inf;\n"
			"      for k=1:3\n"
			"        try\n"
			"          gm = fitgmdist(raw(:),k,'RegularizationValue',1e-6,'Options',statset('MaxIter',500));\n"
			"          BIC = gm.BIC;\n"
			"          if BIC<bestBIC, bestBIC=BIC; bestGM=gm; end\n"
			"        catch, end\n"
			"      end\n"
			"      if ~isempty(bestGM)\n"
			"        xg = linspace(min(bin_centers), max(bin_centers), 400);\n"
			"        yg = pdf(bestGM, xg'); yg = yg(:)';\n"
			"        sfG = max(counts)/max(yg+eps);\n"
			"        plot(xg, yg*sfG, 'm-', 'LineWidth', 1.2);\n"
			"        legend('Histogram','Normal Fit','Mean','KDE','BestParam','GMM','Location','best');\n"
			"      end\n"
			"    end\n"
			"  catch ME\n"
			"    warning('GMM overlay skipped: %s', ME.message);\n"
			"  end\n"
			"end\n\n"
			);
}

/**
 * emit_gmm_overlay_comprehensive
 *
 * Writes an optional GMM overlay section into the comprehensive MATLAB script. It
 * attempts 1–3 component fits with `fitgmdist`, selects the lowest-BIC model, and
 * overlays the scaled PDF for the current field in the loop. All operations are
 * guarded for toolbox availability and sample sufficiency.
 *
 * Works by emitting MATLAB commands that build or load samples from `(bc, cnt)`,
 * fit mixtures with regularization, compare BIC, and overlay the best PDF.
 *
 * @note Expects `fn`, `bc`, and `cnt` to be defined by the surrounding script.
 *
 * @param comp Pointer to an open FILE stream for the comprehensive MATLAB script.
 * @return void.
 */
static void emit_gmm_overlay_comprehensive(FILE *comp)
{
	fprintf(comp,
			"%%%% %% Section: Gaussian Mixture Model (BIC, 1..3 comps)\n"
			"if exist('fitgmdist','file')\n"
			"  try\n"
			"    haveRaw=false; rawCandidates={fullfile(analysisDir,sprintf('%%s_samples.txt',fn)),sprintf('%%s_samples.txt',fn)};\n"
			"    for rk=1:numel(rawCandidates)\n"
			"      if exist(rawCandidates{rk},'file'), try raw=load(rawCandidates{rk}); haveRaw=true; break; end, end\n"
			"    end\n"
			"    if ~haveRaw\n"
			"      raw=[]; for ii=1:numel(bc), raw=[raw; repmat(bc(ii), cnt(ii), 1)]; end\n"
			"    end\n"
			"    if numel(raw)>10\n"
			"      bestGM=[]; bestBIC=Inf;\n"
			"      for k=1:3\n"
			"        try\n"
			"          gm = fitgmdist(raw(:),k,'RegularizationValue',1e-6,'Options',statset('MaxIter',500));\n"
			"          BIC = gm.BIC;\n"
			"          if BIC<bestBIC, bestBIC=BIC; bestGM=gm; end\n"
			"        catch, end\n"
			"      end\n"
			"      if ~isempty(bestGM)\n"
			"        xg = linspace(min(bc), max(bc), 400);\n"
			"        yg = pdf(bestGM, xg'); yg = yg(:)';\n"
			"        sfG = max(cnt)/max(yg+eps);\n"
			"        plot(xg, yg*sfG, 'm-', 'LineWidth', 1.2);\n"
			"        legend('Histogram','Normal Fit','Mean','KDE','BestParam','GMM','Location','best');\n"
			"      end\n"
			"    end\n"
			"  catch ME\n"
			"    warning('GMM overlay skipped: %s', ME.message);\n"
			"  end\n"
			"end\n\n"
			);
}







/**
 * emit_boxcox_scan_individual
 *
 * Writes a Box-Cox transformation scan section for a per-field MATLAB script. When
 * samples are positive, it sweeps λ ∈ [-2, 2] in steps of 0.1, evaluates a normality
 * statistic (Anderson–Darling if available, else KS, else variance proxy), reports the
 * best λ to the console, and leaves plots unchanged.
 *
 * Works by emitting MATLAB code that loads or synthesizes samples, applies the
 * Box-Cox transform across λ, scores normality with available tests, tracks the best
 * score, and prints the selected λ with the field name.
 *
 * @note Requires `analysisDir`, `field`, `bin_centers`, and `counts` to be present.
 *
 * @param mf Pointer to an open FILE stream for the per-field MATLAB script.
 * @return void.
 */
static void emit_boxcox_scan_individual(FILE *mf)
{
	fprintf(mf,
			"%%%% %% Section: Box-Cox scan (report best lambda)\n"
			"bestLam=NaN; bestScore=Inf; haveRaw=false;\n"
			"rawCandidates={fullfile(analysisDir,sprintf('%%s_samples.txt',field)),sprintf('%%s_samples.txt',field)};\n"
			"for rk=1:numel(rawCandidates)\n"
			"  if exist(rawCandidates{rk},'file'), try raw=load(rawCandidates{rk}); haveRaw=true; break; end, end\n"
			"end\n"
			"if ~haveRaw\n"
			"  raw=[]; for ii=1:numel(bin_centers), raw=[raw; repmat(bin_centers(ii), counts(ii), 1)]; end\n"
			"end\n"
			"raw = raw(:);\n"
			"if ~isempty(raw) && all(raw>0)\n"
			"  lam= -2:0.1:2;\n"
			"  for L=lam\n"
			"    y = (L==0) .* log(raw) + (L~=0) .* ((raw.^L - 1)./L);\n"
			"    try\n"
			"      if exist('adtest','file'), [~,p,stat]=adtest((y-mean(y))/std(y)); score=stat; else\n"
			"        if exist('kstest','file'), [~,p]=kstest((y-mean(y))/std(y)); score=1-p; else, score=var(y); end\n"
			"      end\n"
			"      if score<bestScore, bestScore=score; bestLam=L; end\n"
			"    catch, end\n"
			"  end\n"
			"  if ~isnan(bestLam)\n"
			"    disp(sprintf('Box-Cox best lambda for %%s: %%g', field, bestLam));\n"
			"  end\n"
			"end\n\n"
			);
}



/**
 * emit_boxcox_scan_comprehensive
 *
 * Writes a Box-Cox transformation scan section into the comprehensive MATLAB script.
 * For each field in the loop, when data are strictly positive, it sweeps λ, evaluates
 * normality using AD or KS if available, selects the best λ, and prints it to the
 * console. No figures are modified by this step.
 *
 * Works by emitting MATLAB commands that construct/obtain samples from `(bc, cnt)`,
 * transform, evaluate, and report the best λ.
 *
 * @note Expects `fn`, `bc`, and `cnt` to be defined by the caller.
 *
 * @param comp Pointer to an open FILE stream for the comprehensive MATLAB script.
 * @return void.
 */
static void emit_boxcox_scan_comprehensive(FILE *comp)
{
	fprintf(comp,
			"%%%% %% Section: Box-Cox scan (report best lambda)\n"
			"bestLam=NaN; bestScore=Inf; haveRaw=false;\n"
			"rawCandidates={fullfile(analysisDir,sprintf('%%s_samples.txt',fn)),sprintf('%%s_samples.txt',fn)};\n"
			"for rk=1:numel(rawCandidates)\n"
			"  if exist(rawCandidates{rk},'file'), try raw=load(rawCandidates{rk}); haveRaw=true; break; end, end\n"
			"end\n"
			"if ~haveRaw\n"
			"  raw=[]; for ii=1:numel(bc), raw=[raw; repmat(bc(ii), cnt(ii), 1)]; end\n"
			"end\n"
			"raw = raw(:);\n"
			"if ~isempty(raw) && all(raw>0)\n"
			"  lam= -2:0.1:2;\n"
			"  for L=lam\n"
			"    y = (L==0) .* log(raw) + (L~=0) .* ((raw.^L - 1)./L);\n"
			"    try\n"
			"      if exist('adtest','file'), [~,p,stat]=adtest((y-mean(y))/std(y)); score=stat; else\n"
			"        if exist('kstest','file'), [~,p]=kstest((y-mean(y))/std(y)); score=1-p; else, score=var(y); end\n"
			"      end\n"
			"      if score<bestScore, bestScore=score; bestLam=L; end\n"
			"    catch, end\n"
			"  end\n"
			"  if ~isnan(bestLam)\n"
			"    disp(sprintf('Box-Cox best lambda for %%s: %%g', fn, bestLam));\n"
			"  end\n"
			"end\n\n"
			);
}





/**
 * emit_gof_bootstrap_robust_individual
 *
 * Writes a goodness-of-fit / bootstrap / robust-markers section into a per-field
 * MATLAB script. It computes normality p-values (KS/AD if available), bootstrap CIs
 * for mean/std (when enough samples), and overlays robust indicators (median line and
 * IQR band) on the current histogram axes.
 *
 * Works by emitting MATLAB code that loads/synthesizes raw samples, normalizes them,
 * runs available tests, performs bootstrap resampling to estimate CIs, and plots
 * robust summaries with patch/line overlays; any errors are caught and warned.
 *
 * @note Expects `analysisDir`, `field`, `bin_centers`, `counts`, `meanVal`, and
 *       `stdVal` to be defined in the host script.
 *
 * @param mf Pointer to an open FILE stream for the per-field MATLAB script.
 * @return void.
 */
static void emit_gof_bootstrap_robust_individual(FILE *mf)
{
	fprintf(mf,
			"%%%% %% Section: GoF + Bootstrap + Robust\n"
			"%% Normality tests (if available)\n"
			"try\n"
			"  haveRaw=false; rawCandidates={fullfile(analysisDir,sprintf('%%s_samples.txt',field)),sprintf('%%s_samples.txt',field)};\n"
			"  for rk=1:numel(rawCandidates)\n"
			"    if exist(rawCandidates{rk},'file'), try raw=load(rawCandidates{rk}); haveRaw=true; break; end, end\n"
			"  end\n"
			"  if ~haveRaw\n"
			"    raw=[]; for ii=1:numel(bin_centers), raw=[raw; repmat(bin_centers(ii), counts(ii), 1)]; end\n"
			"  end\n"
			"  raw=raw(:);\n"
			"  if ~isempty(raw)\n"
			"    Z = (raw-meanVal)/stdVal;\n"
			"    pKS=NaN; pAD=NaN;\n"
			"    if exist('kstest','file'), [~,pKS]=kstest(Z); end\n"
			"    if exist('adtest','file'), [~,pAD]=adtest(Z); end\n"
			"    disp(sprintf('GoF (normality)  KS p=%%g   AD p=%%g', pKS, pAD));\n"
			"    %% Bootstrap CIs\n"
			"    B=1000; if numel(raw)>10\n"
			"      M=zeros(B,1); S=zeros(B,1);\n"
			"      for b=1:B\n"
			"        idx = randi(numel(raw), numel(raw), 1); rb = raw(idx);\n"
			"        M(b)=mean(rb); S(b)=std(rb);\n"
			"      end\n"
			"      ciM=quantile(M,[0.025 0.975]); ciS=quantile(S,[0.025 0.975]);\n"
			"      disp(sprintf('Bootstrap  mean CI [%%g, %%g],  std CI [%%g, %%g]', ciM(1),ciM(2),ciS(1),ciS(2)));\n"
			"    end\n"
			"    %% Robust markers\n"
			"    med = median(raw); q = quantile(raw,[0.25 0.75]);\n"
			"    yl = ylim; patch([q(1) q(2) q(2) q(1)],[0 0 yl(2) yl(2)], [0.9 0.9 1.0], 'FaceAlpha',0.15, 'EdgeColor','none');\n"
			"    plot([med med],[0 yl(2)],'b-.','LineWidth',1.0);\n"
			"  end\n"
			"catch ME\n"
			"  warning('GoF/Bootstrap skipped: %s',ME.message);\n"
			"end\n\n"
			);
}

/**
 * emit_gof_bootstrap_robust_comprehensive
 *
 * Writes a goodness-of-fit / bootstrap / robust-markers section into the
 * comprehensive MATLAB script. For the current field, it computes KS/AD p-values when
 * available, bootstraps mean/std confidence intervals, and overlays robust markers
 * (median and IQR band) on the active histogram axes.
 *
 * Works by emitting MATLAB commands that assemble samples, run tests, perform
 * bootstrap, and draw robust overlays. Failures are handled via try/catch with a
 * warning message.
 *
 * @note Expects `fn`, `bc`, `cnt`, `mu`, and `sig` to be defined by the caller.
 *
 * @param comp Pointer to an open FILE stream for the comprehensive MATLAB script.
 * @return void.
 */
static void emit_gof_bootstrap_robust_comprehensive(FILE *comp)
{
	fprintf(comp,
			"%%%% %% Section: GoF + Bootstrap + Robust\n"
			"try\n"
			"  haveRaw=false; rawCandidates={fullfile(analysisDir,sprintf('%%s_samples.txt',fn)),sprintf('%%s_samples.txt',fn)};\n"
			"  for rk=1:numel(rawCandidates)\n"
			"    if exist(rawCandidates{rk},'file'), try raw=load(rawCandidates{rk}); haveRaw=true; break; end, end\n"
			"  end\n"
			"  if ~haveRaw\n"
			"    raw=[]; for ii=1:numel(bc), raw=[raw; repmat(bc(ii), cnt(ii), 1)]; end\n"
			"  end\n"
			"  raw=raw(:);\n"
			"  if ~isempty(raw)\n"
			"    Z = (raw-mu)/sig;\n"
			"    pKS=NaN; pAD=NaN;\n"
			"    if exist('kstest','file'), [~,pKS]=kstest(Z); end\n"
			"    if exist('adtest','file'), [~,pAD]=adtest(Z); end\n"
			"    disp(sprintf('GoF (normality)  KS p=%%g   AD p=%%g', pKS, pAD));\n"
			"    B=1000; if numel(raw)>10\n"
			"      M=zeros(B,1); S=zeros(B,1);\n"
			"      for b=1:B\n"
			"        idx = randi(numel(raw), numel(raw), 1); rb = raw(idx);\n"
			"        M(b)=mean(rb); S(b)=std(rb);\n"
			"      end\n"
			"      ciM=quantile(M,[0.025 0.975]); ciS=quantile(S,[0.025 0.975]);\n"
			"      disp(sprintf('Bootstrap  mean CI [%%g, %%g],  std CI [%%g, %%g]', ciM(1),ciM(2),ciS(1),ciS(2)));\n"
			"    end\n"
			"    med=median(raw); q=quantile(raw,[0.25 0.75]); yl=ylim;\n"
			"    patch([q(1) q(2) q(2) q(1)],[0 0 yl(2) yl(2)],[0.9 0.9 1.0],'FaceAlpha',0.15,'EdgeColor','none');\n"
			"    plot([med med],[0 yl(2)],'b-.','LineWidth',1.0);\n"
			"  end\n"
			"catch ME\n"
			"  warning('GoF/Bootstrap skipped: %s',ME.message);\n"
			"end\n\n"
			);
}













/* ========================================================================
 * Internal helpers (static)
 * ===================================================================== */
/**
 * resolve_scripts_dir
 *
 * Computes the directory where MATLAB scripts will be written and ensures it exists
 * when a subdirectory is requested. When `scriptsSubdir` is empty or NULL, the output
 * directory is simply `analysisDir`. The resulting path is written into `dst`.
 *
 * Works by formatting either `analysisDir/scriptsSubdir` (and creating it) or the
 * bare `analysisDir` into the destination buffer. The caller then uses `dst` for
 * subsequent file creations.
 *
 * @note If `create_directory` fails, later file operations will report errors; this
 *       function does not abort on directory creation failure.
 *
 * @param dst Buffer to receive the resolved directory path (NUL-terminated).
 * @param dstSize Size of `dst` in bytes.
 * @param analysisDir Base directory that already exists from previous pipeline steps.
 * @param scriptsSubdir Optional subdirectory name under `analysisDir` for outputs.
 * @return void.
 */
static void resolve_scripts_dir(char *dst, size_t dstSize,
								const char *analysisDir,
								const char *scriptsSubdir)
{
	/// Decide target: either analysisDir/<scriptsSubdir> or analysisDir itself
	if (scriptsSubdir && scriptsSubdir[0] != '\0') {
		snprintf(dst, dstSize, "%s/%s", analysisDir, scriptsSubdir);  // Compose path
		create_directory(dst, "");                                     // Ensure it exists
	} else {
		snprintf(dst, dstSize, "%s", analysisDir);                     // Use base directly
																	   // analysisDir is assumed to exist already.
	}
}



/**
 * emit_individual_script_for_field
 *
 * Generates a self-contained per-field MATLAB script `<field>_plot.m` that robustly
 * discovers histogram/statistics files, plots the histogram, overlays a scaled
 * Normal fit with a mean marker, prints a concise console summary, and (optionally)
 * appends advanced modeling/KDE/GMM/Box-Cox/GoF sections. A PNG of the figure is saved
 * next to the script at MATLAB runtime.
 *
 * Works by opening `<outDir>/<fieldName>_plot.m` and writing MATLAB code that:
 * (1) resolves `histogram` and `stats` via multiple candidate paths, (2) builds the
 * baseline plot and legend, (3) emits optional modeling snippets, and (4) saves a PNG.
 * The function reports file-creation errors and returns without emitting on failure.
 *
 * @note The emitted MATLAB script is resilient to missing files/toolboxes and uses
 *       try/catch and `exist(...)` guards accordingly.
 *
 * @param outDir Directory where the MATLAB script is created.
 * @param analysisDir Preferred directory to locate input files at runtime.
 * @param fieldName Logical field name used to resolve file names and labels.
 * @param indexForFigSuffix Integer suffix to make MATLAB figure handles unique.
 * @return void.
 */
static void emit_individual_script_for_field(const char *outDir,
											 const char *analysisDir,
											 const char *fieldName,
											 int indexForFigSuffix)
{
	char scriptPath[4096];
	snprintf(scriptPath, sizeof(scriptPath), "%s/%s_plot.m", outDir, fieldName);
	
	FILE *mf = fopen(scriptPath, "w");
	if (!mf) {
		fprintf(stderr, "Error creating MATLAB plot script: %s\n", scriptPath);
		return;
	}
	
	/* Header + core inputs (keep original intent) */
	fprintf(mf, "%%%% Auto-generated: per-field modeling for \"%s\"\n", fieldName);
	fprintf(mf, "field       = '%s';\n", fieldName);
	fprintf(mf, "analysisDir = '%s';\n", analysisDir ? analysisDir : "");
	/* Resolve the directory this script lives in (so saves land next to the .m file) */
	fprintf(mf, "scriptDir = fileparts(mfilename('fullpath'));\n\n");
	
	/* --- Histogram discovery (compatible superset) ---
	 Prefer the original file in analysisDir; fall back to scriptDir or CWD. */
	fprintf(mf, "histCandidates = { ...\n");
	fprintf(mf, "  fullfile(analysisDir, sprintf('%%s_histogram.txt', field)), ...\n");
	fprintf(mf, "  fullfile(scriptDir,   sprintf('%%s_histogram.txt', field)), ...\n");
	fprintf(mf, "  sprintf('%%s_histogram.txt', field) ...\n");
	fprintf(mf, "};\n");
	fprintf(mf, "histLoaded = false;\n");
	fprintf(mf, "for k = 1:numel(histCandidates)\n");
	fprintf(mf, "    p = histCandidates{k};\n");
	fprintf(mf, "    if exist(p,'file')\n");
	fprintf(mf, "        try\n");
	fprintf(mf, "            hist_data = load(p);\n");
	fprintf(mf, "            if size(hist_data,2) >= 3\n");
	fprintf(mf, "                bin_centers = (hist_data(:,1) + hist_data(:,2))/2;\n");
	fprintf(mf, "                counts      = hist_data(:,3);\n");
	fprintf(mf, "                histLoaded  = true; break;\n");
	fprintf(mf, "            end\n");
	fprintf(mf, "        catch, end\n");
	fprintf(mf, "    end\n");
	fprintf(mf, "end\n");
	fprintf(mf, "if ~histLoaded, error('Histogram not found for %%s', field); end\n\n");
	
	/* Figure + histogram (preserve your figure-handle style) */
	fprintf(mf, "f%d = figure('Name', '%s Histogram');\n", indexForFigSuffix, fieldName);
	fprintf(mf, "bar(bin_centers, counts, 1.0, 'FaceColor',[0.7 0.7 0.7]); grid on; hold on;\n");
	fprintf(mf, "title('%s Histogram'); xlabel('%s Values'); ylabel('Frequency');\n\n", fieldName, fieldName);

	/* --- Stats discovery (compatible superset) ---
	 Try *_stats.txt, *_full_analysis.txt, *_analysis.txt in analysisDir, then scriptDir, then CWD. */
	fprintf(mf, "statNames = { '_stats.txt', '_full_analysis.txt', '_analysis.txt' };\n");
	fprintf(mf, "statRoots = { analysisDir, scriptDir, '' };\n");
	fprintf(mf, "statsTxt = '';\n");
	fprintf(mf, "for r = 1:numel(statRoots)\n");
	fprintf(mf, "  for n = 1:numel(statNames)\n");
	fprintf(mf, "    if isempty(statRoots{r}), p = sprintf('%%s%%s', field, statNames{n});\n");
	fprintf(mf, "    else, p = fullfile(statRoots{r}, sprintf('%%s%%s', field, statNames{n})); end\n");
	fprintf(mf, "    if exist(p,'file'), statsTxt = fileread(p); break; end\n");
	fprintf(mf, "  end\n");
	fprintf(mf, "  if ~isempty(statsTxt), break; end\n");
	fprintf(mf, "end\n");
	
	/* Parse mean/std (your originals), plus optional skew & AD (safe extras) */
	fprintf(mf, "meanVal = 0; stdVal = 1; skewVal = NaN; adVal = NaN;\n");
	fprintf(mf, "if ~isempty(statsTxt)\n");
	fprintf(mf, "  tokens_mean = regexp(statsTxt, 'Mean:\\s*([-+0-9eE\\.]+)', 'tokens');\n");
	fprintf(mf, "  if ~isempty(tokens_mean), meanVal = str2double(tokens_mean{1}{1}); end\n");
	fprintf(mf, "  tstd = regexp(statsTxt, '(?:Std Dev|Standard Deviation):\\s*([-+0-9eE\\.]+)', 'tokens');\n");
	fprintf(mf, "  if ~isempty(tstd), stdVal = str2double(tstd{1}{1}); end\n");
	fprintf(mf, "  tsk  = regexp(statsTxt, 'Skewness:\\s*([-+0-9eE\\.]+)', 'tokens');\n");
	fprintf(mf, "  if ~isempty(tsk),  skewVal = str2double(tsk{1}{1}); end\n");
	fprintf(mf, "  tad  = regexp(statsTxt, 'Anderson-Darling A\\^2:\\s*([-+0-9eE\\.]+)', 'tokens');\n");
	fprintf(mf, "  if ~isempty(tad),  adVal   = str2double(tad{1}{1}); end\n");
	fprintf(mf, "end\n\n");
	
	/* Normal overlay (keep your scaling using eps) + mean marker + legend */
	fprintf(mf, "x  = linspace(min(bin_centers), max(bin_centers), 200);\n");
	fprintf(mf, "y  = (1./(stdVal.*sqrt(2*pi))) .* exp(-0.5*((x-meanVal)./stdVal).^2);\n");
	fprintf(mf, "sf = max(counts) / max(y + eps);\n");
	fprintf(mf, "plot(x, y*sf, 'r-', 'LineWidth', 1.5);\n");
	fprintf(mf, "plot([meanVal meanVal], [0 max(counts)], 'k--', 'LineWidth', 1.2);\n");
	fprintf(mf, "legend('Histogram','Normal Fit','Mean');\n\n");
	
	/* Console summary (safe optional extras) */
	fprintf(mf, "disp('======================');\n");
	fprintf(mf, "disp(['Field: ', field]);\n");
	fprintf(mf, "disp(['Mean = ', num2str(meanVal)]);\n");
	fprintf(mf, "disp(['Std  = ', num2str(stdVal)]);\n");
	fprintf(mf, "if ~isnan(skewVal), disp(['Skew = ', num2str(skewVal)]); end\n");
	fprintf(mf, "if ~isnan(adVal),   disp(['Anderson-Darling = ', num2str(adVal)]); end\n");
	fprintf(mf, "disp('======================');\n\n");
	
	
	emit_adv_modeling_snippet_individual(mf);
	emit_kde_overlay_individual(mf);
	emit_gmm_overlay_individual(mf);
	emit_boxcox_scan_individual(mf);
	emit_gof_bootstrap_robust_individual(mf);
	
	
	/* Save PNG next to the .m script (non-intrusive improvement) */
	fprintf(mf, "saveas(gcf, fullfile(scriptDir, sprintf('%%s_plot.png', field)));\n");
	
	fclose(mf);
}



/**
 * emit_comprehensive_script
 *
 * Generates a single MATLAB script `comprehensive_plots.m` that iterates over all
 * fields, robustly discovers each field’s histogram and stats, produces the baseline
 * histogram + scaled Normal overlay + mean marker, prints a summary to the console,
 * appends optional modeling/KDE/GMM/Box-Cox/GoF sections, and saves a PNG per field.
 *
 * Works by creating `<outDir>/comprehensive_plots.m` and emitting MATLAB code that
 * sets paths, defines the `fields` cell array from `fieldNames`, loops through each
 * field with guarded discovery, plotting, optional overlays, and per-field saves.
 * File-creation errors are reported and abort emission.
 *
 * @param outDir Directory where the comprehensive MATLAB script is written.
 * @param analysisDir Preferred location for input discovery at runtime.
 * @param fieldNames Array of field name strings to enumerate in the script.
 * @param numFields Number of entries in `fieldNames`.
 * @return void.
 */
static void emit_comprehensive_script(const char *outDir,
									  const char *analysisDir,
									  char **fieldNames,
									  int numFields)
{
	char compPath[4096];
	snprintf(compPath, sizeof(compPath), "%s/comprehensive_plots.m", outDir);
	
	FILE *comp = fopen(compPath, "w");
	if (!comp) {
		fprintf(stderr, "Error creating comprehensive_plots.m\n");
		return;
	}
	
	fprintf(comp, "%%%% Auto-generated: comprehensive modeling over all fields\n");
	fprintf(comp, "close all; clear; clc;\n\n");
	
	/* Script-local paths */
	fprintf(comp, "analysisDir = '%s';\n", analysisDir ? analysisDir : "");
	fprintf(comp, "scriptDir   = fileparts(mfilename('fullpath'));\n\n");
	
	/* Emit fields cell array */
	fprintf(comp, "fields = {");
	for (int i = 0; i < numFields; i++)
	{
		if (fieldNames[i] && fieldNames[i][0] != '\\0') {
			fprintf(comp, "'%s'%s", fieldNames[i], (i < numFields - 1) ? ", " : "");
		}
	}
	fprintf(comp, "};\n\n");
	
	/* Unified discovery & plotting per field (non-fatal on missing data) */
	fprintf(comp, "for i = 1:numel(fields)\n");
	fprintf(comp, "  fn = fields{i};\n");
	fprintf(comp, "  try\n");
	fprintf(comp, "    %% --- Histogram discovery (analysisDir → scriptDir → CWD) ---\n");
	fprintf(comp, "    histCandidates = {\n");
	fprintf(comp, "      fullfile(analysisDir, sprintf('%%s_histogram.txt', fn)),\n");
	fprintf(comp, "      fullfile(scriptDir,   sprintf('%%s_histogram.txt', fn)),\n");
	fprintf(comp, "      sprintf('%%s_histogram.txt', fn)\n");
	fprintf(comp, "    };\n");
	fprintf(comp, "    histLoaded = false;\n");
	fprintf(comp, "    for k = 1:numel(histCandidates)\n");
	fprintf(comp, "      p = histCandidates{k};\n");
	fprintf(comp, "      if exist(p,'file')\n");
	fprintf(comp, "        try\n");
	fprintf(comp, "          H = load(p);\n");
	fprintf(comp, "          if size(H,2) >= 3\n");
	fprintf(comp, "            bc  = (H(:,1)+H(:,2))/2;\n");
	fprintf(comp, "            cnt = H(:,3);\n");
	fprintf(comp, "            histLoaded = true; break;\n");
	fprintf(comp, "          end\n");
	fprintf(comp, "        catch, end\n");
	fprintf(comp, "      end\n");
	fprintf(comp, "    end\n");
	fprintf(comp, "    if ~histLoaded\n");
	fprintf(comp, "      warning('Histogram not found for %%s. Skipping.', fn);\n");
	fprintf(comp, "      continue;\n");
	fprintf(comp, "    end\n\n");
	
	fprintf(comp, "    %% --- Stats discovery: *_stats|_full_analysis|_analysis in analysisDir/scriptDir/CWD ---\n");
	fprintf(comp, "    statNames = { '_stats.txt', '_full_analysis.txt', '_analysis.txt' };\n");
	fprintf(comp, "    statRoots = { analysisDir, scriptDir, '' };\n");
	fprintf(comp, "    statsTxt = '';\n");
	fprintf(comp, "    for r = 1:numel(statRoots)\n");
	fprintf(comp, "      for n = 1:numel(statNames)\n");
	fprintf(comp, "        if isempty(statRoots{r})\n");
	fprintf(comp, "          p = sprintf('%%s%%s', fn, statNames{n});\n");
	fprintf(comp, "        else\n");
	fprintf(comp, "          p = fullfile(statRoots{r}, sprintf('%%s%%s', fn, statNames{n}));\n");
	fprintf(comp, "        end\n");
	fprintf(comp, "        if exist(p,'file'), statsTxt = fileread(p); break; end\n");
	fprintf(comp, "      end\n");
	fprintf(comp, "      if ~isempty(statsTxt), break; end\n");
	fprintf(comp, "    end\n");
	fprintf(comp, "    mu = 0; sig = 1; skewVal = NaN; adVal = NaN;\n");
	fprintf(comp, "    if ~isempty(statsTxt)\n");
	fprintf(comp, "      tm = regexp(statsTxt, 'Mean:\\s*([-+0-9eE\\.]+)', 'tokens');\n");
	fprintf(comp, "      if ~isempty(tm), mu = str2double(tm{1}{1}); end\n");
	fprintf(comp, "      ts = regexp(statsTxt, '(?:Std Dev|Standard Deviation):\\s*([-+0-9eE\\.]+)', 'tokens');\n");
	fprintf(comp, "      if ~isempty(ts), sig = str2double(ts{1}{1}); end\n");
	fprintf(comp, "      tsk = regexp(statsTxt, 'Skewness:\\s*([-+0-9eE\\.]+)', 'tokens');\n");
	fprintf(comp, "      if ~isempty(tsk), skewVal = str2double(tsk{1}{1}); end\n");
	fprintf(comp, "      tad = regexp(statsTxt, 'Anderson-Darling A\\^2:\\s*([-+0-9eE\\.]+)', 'tokens');\n");
	fprintf(comp, "      if ~isempty(tad), adVal = str2double(tad{1}{1}); end\n");
	fprintf(comp, "    end\n");
	fprintf(comp, "    if ~(isfinite(sig) && sig > 0), sig = 1; end\n\n");
	
	fprintf(comp, "    %% --- Plot: histogram + scaled normal + mean marker ---\n");
	fprintf(comp, "    f = figure('Name',[fn ' Histogram']);\n");
	fprintf(comp, "    bar(bc, cnt, 1.0, 'FaceColor',[0.7 0.7 0.7]); grid on; hold on;\n");
	fprintf(comp, "    title([fn ' Histogram']); xlabel([fn ' Values']); ylabel('Frequency');\n");
	fprintf(comp, "    x  = linspace(min(bc), max(bc), 200);\n");
	fprintf(comp, "    y  = (1./(sig.*sqrt(2*pi))) .* exp(-0.5*((x-mu)./sig).^2);\n");
	fprintf(comp, "    sf = max(cnt) / max(y + eps);\n");
	fprintf(comp, "    plot(x, y*sf, 'r-', 'LineWidth', 1.5);\n");
	fprintf(comp, "    plot([mu mu], [0 max(cnt)], 'k--', 'LineWidth', 1.2);\n");
	fprintf(comp, "    legend('Histogram','Normal Fit','Mean','Location','best');\n\n");
	
	fprintf(comp, "    %% --- Console summary ---\n");
	fprintf(comp, "    disp('======================');\n");
	fprintf(comp, "    disp(['Field: ', fn]);\n");
	fprintf(comp, "    disp(['Mean = ', num2str(mu)]);\n");
	fprintf(comp, "    disp(['Std  = ', num2str(sig)]);\n");
	fprintf(comp, "    if ~isnan(skewVal), disp(['Skew = ', num2str(skewVal)]); end\n");
	fprintf(comp, "    if ~isnan(adVal),   disp(['Anderson-Darling = ', num2str(adVal)]); end\n");
	fprintf(comp, "    disp('======================');\n\n");
	
	emit_adv_modeling_snippet_comprehensive(comp);
	emit_kde_overlay_comprehensive(comp);
	emit_gmm_overlay_comprehensive(comp);
	emit_boxcox_scan_comprehensive(comp);
	emit_gof_bootstrap_robust_comprehensive(comp);
	
	
	fprintf(comp, "    %% --- Save PNG beside this script ---\n");
	fprintf(comp, "    saveas(f, fullfile(scriptDir, sprintf('%%s_comprehensive.png', fn)));\n");
	fprintf(comp, "  catch ME\n");
	fprintf(comp, "    warning('Error processing %%s: %%s', fn, ME.message);\n");
	fprintf(comp, "  end\n");
	fprintf(comp, "end\n");
	
	fclose(comp);
}









/**
 * generate_matlab_scripts_unified
 *
 * Public entry that emits MATLAB plotting/modeling scripts in one pass. Depending on
 * `flavor`, it generates per-field scripts `<field>_plot.m`, the aggregate
 * `comprehensive_plots.m`, or both, placing them under `analysisDir` or an optional
 * `scriptsSubdir`.
 *
 * Works by validating inputs, resolving the output directory, iterating fields to
 * emit individual scripts when requested, and emitting the comprehensive script when
 * requested. Failures to create specific files are reported per case without aborting
 * the entire generation.
 *
 * @param analysisDir Base directory used for input discovery within the MATLAB code.
 * @param fieldNames Array of field name strings to target.
 * @param numFields Number of field names.
 * @param flavor Bitmask selecting INDIVIDUAL and/or COMPREHENSIVE outputs.
 * @param scriptsSubdir Optional subfolder under `analysisDir` for script outputs.
 * @return void.
 */
void generate_matlab_scripts_unified(const char *analysisDir,
									 char **fieldNames,
									 int numFields,
									 MatlabScriptFlavor flavor,
									 const char *scriptsSubdir)
{
	if (!analysisDir || !fieldNames || numFields <= 0) {
		fprintf(stderr, "generate_matlab_scripts_unified: invalid arguments.\n");
		return;
	}
	
	char outDir[4096];
	resolve_scripts_dir(outDir, sizeof(outDir), analysisDir, scriptsSubdir);
	
	if (flavor & MATLAB_SCRIPTS_INDIVIDUAL) {
		for (int i = 0; i < numFields; ++i) {
			if (fieldNames[i] && fieldNames[i][0] != '\0') {
				emit_individual_script_for_field(outDir, analysisDir, fieldNames[i], i);
			}
		}
	}
	
	if (flavor & MATLAB_SCRIPTS_COMPREHENSIVE) {
		emit_comprehensive_script(outDir, analysisDir, fieldNames, numFields);
	}
}





































































/**
 * generate_matlab_model_scripts_unified
 *
 * Unified emitter for standardized per-field modeling scripts `<field>_model.m`
 * written to `<analysisDir>/MATLAB_Scripts`. Each MATLAB script robustly discovers
 * the histogram and stats files, plots the histogram with a scaled Normal overlay,
 * annotates/prints a concise summary, and saves `<field>_modeling.png` at runtime.
 *
 * Works by ensuring `<analysisDir>/MATLAB_Scripts` exists, then, for each field,
 * opening/creating `<field>_model.m` (append or overwrite as requested) and writing
 * MATLAB code that: (1) discovers inputs from multiple candidate paths, (2) builds the
 * baseline plot and legend, (3) prints a summary, and (4) saves the PNG. Errors are
 * logged per field without aborting the batch.
 *
 * @note The emitted MATLAB is guarded with `exist(...)` checks and try/catch blocks
 *       to remain operational in GNU Octave or MATLAB without specific toolboxes.
 *
 * @param analysisDir Base directory (also the parent for the MATLAB_Scripts folder).
 * @param fieldNames Array of field name strings to generate scripts for.
 * @param numFields Number of entries in `fieldNames`.
 * @param overwrite If 0, append to existing file; if nonzero, truncate and rewrite.
 * @return void.
 */
void generate_matlab_model_scripts_unified(const char *analysisDir,
										   char **fieldNames,
										   int numFields,
										   int overwrite)   // 0 = append, nonzero = overwrite
{
	// Choose target parent dir (prefer analysisDir if provided; else resultsDir)
	const char *parent = analysisDir;
	
	
	// scripts output folder: <parent>/MATLAB_Scripts
	char scriptsDir[4096];
	snprintf(scriptsDir, sizeof(scriptsDir), "%s/%s", parent, "MATLAB_Scripts");
	create_directory(scriptsDir, "");
	
	for (int i = 0; i < numFields; i++) {
		const char *field = fieldNames[i];
		if (!field || !*field) continue;
		
		char scriptPath[4096];
		snprintf(scriptPath, sizeof(scriptPath), "%s/%s_model.m", scriptsDir, field);
		
		FILE *mf = fopen(scriptPath, overwrite ? "w" : "a+");
		if (!mf) {
			fprintf(stderr, "Error creating MATLAB script: %s\n", scriptPath);
			continue;
		}
		
		// ---------- MATLAB SCRIPT ----------
		fprintf(mf, "%%%% Auto-generated unified modeling script\n");
		fprintf(mf, "field       = '%s';\n", field);
		fprintf(mf, "analysisDir = '%s';\n", analysisDir ? analysisDir : "");
		
		// Locate histogram
		fprintf(mf, "histCandidates = {\n");
		fprintf(mf, "  fullfile(analysisDir, sprintf('%%s_histogram.txt', field)),\n");
		fprintf(mf, "  sprintf('%%s_histogram.txt', field)\n");
		fprintf(mf, "};\n");
		fprintf(mf, "histLoaded = false;\n");
		fprintf(mf, "for k = 1:numel(histCandidates)\n");
		fprintf(mf, "    p = histCandidates{k}; if exist(p,'file')\n");
		fprintf(mf, "        try\n");
		fprintf(mf, "            H = load(p);\n");
		fprintf(mf, "            if size(H,2) >= 3\n");
		fprintf(mf, "                bin_centers = (H(:,1)+H(:,2))/2;\n");
		fprintf(mf, "                counts      = H(:,3);\n");
		fprintf(mf, "                histLoaded  = true; break;\n");
		fprintf(mf, "            end\n");
		fprintf(mf, "        catch, end\n");
		fprintf(mf, "    end\n");
		fprintf(mf, "end\n");
		fprintf(mf, "if ~histLoaded, error('Histogram not found for %%s', field); end\n\n");
		
		// Locate stats/analysis (try all three names across all dirs)
		fprintf(mf, "meanVal = NaN; stdVal = NaN; skewVal = NaN; adVal = NaN;\n");
		fprintf(mf, "statNames = { '_full_analysis.txt', '_analysis.txt', '_stats.txt' };\n");
		fprintf(mf, "statRoots = { analysisDir, '' };\n");
		fprintf(mf, "for r = 1:numel(statRoots)\n");
		fprintf(mf, "  for n = 1:numel(statNames)\n");
		fprintf(mf, "    cand = statRoots{r};\n");
		fprintf(mf, "    if isempty(cand), p = sprintf('%%s%%s', field, statNames{n});\n");
		fprintf(mf, "    else,              p = fullfile(cand, sprintf('%%s%%s', field, statNames{n})); end\n");
		fprintf(mf, "    if exist(p,'file')\n");
		fprintf(mf, "      txt = fileread(p);\n");
		fprintf(mf, "      t = regexp(txt,'Mean:\\s*([-+0-9eE\\.]+)','tokens'); if ~isempty(t) && isnan(meanVal), meanVal=str2double(t{1}{1}); end\n");
		fprintf(mf, "      t = regexp(txt,'Standard Deviation:\\s*([-+0-9eE\\.]+)','tokens'); if ~isempty(t) && isnan(stdVal), stdVal=str2double(t{1}{1}); end\n");
		fprintf(mf, "      t = regexp(txt,'Std Dev:\\s*([-+0-9eE\\.]+)','tokens');        if ~isempty(t) && isnan(stdVal), stdVal=str2double(t{1}{1}); end\n");
		fprintf(mf, "      t = regexp(txt,'Skewness:\\s*([-+0-9eE\\.]+)','tokens');         if ~isempty(t) && isnan(skewVal), skewVal=str2double(t{1}{1}); end\n");
		fprintf(mf, "      t = regexp(txt,'Anderson-Darling A\\^2:\\s*([-+0-9eE\\.]+)','tokens'); if ~isempty(t) && isnan(adVal), adVal=str2double(t{1}{1}); end\n");
		fprintf(mf, "    end\n");
		fprintf(mf, "  end\n");
		fprintf(mf, "end\n");
		fprintf(mf, "if isnan(meanVal), meanVal = 0; end\n");
		fprintf(mf, "if isnan(stdVal) || stdVal<=0, stdVal = 1; end\n\n");
		
		// Plotting (hist + scaled normal + mean line)
		fprintf(mf, "f = figure('Name', sprintf('%%s Histogram', field));\n");
		fprintf(mf, "bar(bin_centers, counts, 1.0, 'FaceColor',[0.7 0.7 0.7]); hold on; grid on;\n");
		fprintf(mf, "title(sprintf('%%s Histogram', field)); xlabel([field ' Values']); ylabel('Frequency');\n");
		fprintf(mf, "x = linspace(min(bin_centers), max(bin_centers), 200);\n");
		fprintf(mf, "y = (1./(stdVal*sqrt(2*pi))) .* exp(-0.5*((x-meanVal)./stdVal).^2);\n");
		fprintf(mf, "max_counts = max(counts); max_y = max(y); scale = (max_y>0)*(max_counts/max_y) + (max_y==0)*1;\n");
		fprintf(mf, "plot(x, y*scale, 'r-', 'LineWidth', 1.5);\n");
		fprintf(mf, "legend('Histogram','Normal Fit','Location','best');\n");
		fprintf(mf, "plot([meanVal meanVal],[0 max_counts],'k--','LineWidth',1.2);\n");
		fprintf(mf, "text(meanVal, max_counts*0.9, sprintf('Mean=%%.4g',meanVal), 'HorizontalAlignment','center','Color','k');\n\n");
		
		// Console summary + save PNG
		fprintf(mf, "disp('======================');\n");
		fprintf(mf, "disp(['Field: ', field]);\n");
		fprintf(mf, "disp(['Mean = ', num2str(meanVal)]);\n");
		fprintf(mf, "disp(['Std  = ', num2str(stdVal)]);\n");
		fprintf(mf, "if ~isnan(skewVal), disp(['Skew = ', num2str(skewVal)]); end\n");
		fprintf(mf, "if ~isnan(adVal),   disp(['Anderson-Darling = ', num2str(adVal)]); end\n");
		fprintf(mf, "disp('======================');\n");
		fprintf(mf, "saveas(f, sprintf('%%s_modeling.png', field));\n");
		// ---------- /MATLAB SCRIPT ----------
		
		fclose(mf);
	}
}






