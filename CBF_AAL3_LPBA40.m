
%% ============================================================
% NATIVE SPACE CBF BATCH EXTRACTION
% SPM25 / MATLAB
%
% 40 M0 + 40 M1 + 42 H = 122 scans
% Exclude patients: 002, 010, 053, 060
%
% LPBA40: 55 ROI
% AAL3: 150 selected ROI
% GM threshold >= 0.30
%
% No Excel or MAT input dependency
% No custom functions
% Original images unchanged
%% ============================================================

clear; clc;

addpath('/Users/mac/spm','-begin');
rehash;
spm('Defaults','fMRI');
spm_jobman('initcfg');

%% 1. Paths

aslDir = '/Users/mac/EXP/eTRE/ASL';
runDir = fullfile(aslDir,'asl_run');

lpba = ['/Users/mac/LPBA40SPM5/maxprob/' ...
    'lpba40.spm5.avg152T1.gm.label.nii'];

aal = '/Users/mac/AAL3/AAL3v1_1mm.nii';

assert(isfolder(runDir),'ASL扫描目录不存在');
assert(isfile(lpba),'LPBA40模板不存在');
assert(isfile(aal),'AAL3模板不存在');

fprintf('\n输入目录检查通过。\n');

%% 2. Automatically find M0, M1 and H directories

D0 = dir(fullfile(runDir,'*_M0'));
D1 = dir(fullfile(runDir,'*_M1'));
DH = dir(fullfile(runDir,'H*'));

names0 = string({D0([D0.isdir]).name})';
names1 = string({D1([D1.isdir]).name})';
namesH = string({DH([DH.isdir]).name})';

m0 = erase(names0,'_M0');
m1 = erase(names1,'_M1');

m0 = unique(m0(strlength(m0)==3 & ...
    ~isnan(str2double(m0))));

m1 = unique(m1(strlength(m1)==3 & ...
    ~isnan(str2double(m1))));

namesH = unique(namesH( ...
    strlength(namesH)==4 & ...
    ~isnan(str2double(extractAfter(namesH,1)))));

assert(numel(m0)==44 && numel(m1)==44, ...
    '没有检测到完整的44个M0和44个M1目录');

%% 3. Exclusions

excludePatients = ["002";"010";"053";"060"];

excludeH = ["H014";"H018";"H021"; ...
            "H032";"H052";"H072"];

patients = intersect(m0,m1);
patients = sort(setdiff(patients,excludePatients));

healthy = sort(setdiff(namesH,excludeH));

assert(numel(patients)==40, ...
    '患者数量不是40对');

assert(numel(healthy)==42, ...
    '正常对照数量不是42');

%% 4. Build 122-scan list

scan = strings(122,1);
group = strings(122,1);
pid = strings(122,1);

for k = 1:40

    scan(2*k-1) = patients(k)+"_M0";
    scan(2*k)   = patients(k)+"_M1";

    group(2*k-1:2*k) = ["M0";"M1"];
    pid(2*k-1:2*k) = patients(k);

end

scan(81:122) = healthy;
group(81:122) = "H";
pid(81:122) = healthy;

assert(numel(unique(scan))==122, ...
    '扫描编号存在重复');

fprintf('\n========== SAMPLE SUMMARY ==========\n');
fprintf('M0: %d\n',nnz(group=="M0"));
fprintf('M1: %d\n',nnz(group=="M1"));
fprintf('H:  %d\n',nnz(group=="H"));
fprintf('Total: %d\n',numel(scan));

%% 5. Precheck all required images

missingFiles = strings(0,1);

for s = 1:122

    root = fullfile(runDir,char(scan(s)));

    files = {
        fullfile(root,'STRUC','T1.nii')
        fullfile(root,'FUNC', ...
            'meanCBF_0_sASLflt_rPASL.nii')
        fullfile(root,'STRUC', ...
            'segmentation_work','c1T1.nii')
        fullfile(root,'STRUC', ...
            'segmentation_work','y_T1.nii')
    };

    for k = 1:numel(files)

        if ~isfile(files{k})
            missingFiles(end+1,1) = string(files{k});
        end

    end
end

if ~isempty(missingFiles)

    fprintf('\n以下文件缺失：\n');
    disp(missingFiles);

    error('输入文件检查失败，尚未开始批量处理。');

end

fprintf('\n122次扫描的T1、CBF、GM、形变场全部存在。\n');

%% 6. ROI IDs

LPBA_ID = [21:34 41:50 61:68 81:92 ...
           101 102 121 122 161:166 181]';

AAL_ID = [1:34 37:80 83:124 ...
          126:132 135:140 142:158]';

assert(numel(LPBA_ID)==55);
assert(numel(AAL_ID)==150);

%% 7. Output folder

stamp = datestr(now,'yyyymmdd_HHMMSS');

outDir = fullfile(aslDir, ...
    ['NativeCBF_122Scans_' stamp]);

mkdir(outDir);

fprintf('\n输出目录：\n%s\n',outDir);

%% 8. Allocate results

CBF_LPBA = nan(122,55);
CBF_AAL = nan(122,150);

GM_N_LPBA = zeros(122,55);
GM_N_AAL = zeros(122,150);

Finite_N_LPBA = zeros(122,55);
Finite_N_AAL = zeros(122,150);

success = false(122,1);
NeedsReview = false(122,1);
errorMsg = strings(122,1);

Translation_mm = nan(122,3);
Rotation_deg = nan(122,3);

%% 9. Coregistration settings

flags.cost_fun = 'nmi';
flags.sep = [4 2];
flags.tol = [.02 .02 .02 .001 .001 .001];
flags.fwhm = [7 7];
flags.graphics = false;

%% 10. Main batch

for s = 1:122

    sid = char(scan(s));

    fprintf('\n========== [%d/122] %s ==========\n',s,sid);

    try

        root = fullfile(runDir,sid);

        t1 = fullfile(root,'STRUC','T1.nii');

        rawCBF = fullfile(root,'FUNC', ...
            'meanCBF_0_sASLflt_rPASL.nii');

        gmFile = fullfile(root,'STRUC', ...
            'segmentation_work','c1T1.nii');

        yFile = fullfile(root,'STRUC', ...
            'segmentation_work','y_T1.nii');

        subDir = fullfile(outDir,sid);
        mkdir(subDir);

        %% A. Copy and coregister CBF to T1

        cbfCopy = fullfile(subDir, ...
            'Native_Coreg_CBF.nii');

        [ok,msg] = copyfile(rawCBF,cbfCopy);
        assert(ok,'%s',msg);

        Vt = spm_vol(t1);
        Vc = spm_vol(cbfCopy);

        assert(numel(Vt)==1 && numel(Vc)==1, ...
            'T1或CBF不是3D图像');

        oldMat = spm_get_space(cbfCopy);

        x = spm_coreg(Vt,Vc,flags);

        assert(all(isfinite(x(1:6))), ...
            '配准参数存在非有限值');

        spm_get_space(cbfCopy, ...
            spm_matrix(x)\oldMat);

        Translation_mm(s,:) = x(1:3);
        Rotation_deg(s,:) = x(4:6)*180/pi;

        NeedsReview(s) = ...
            any(abs(x(1:3))>30) || ...
            any(abs(x(4:6)*180/pi)>20);

        %% B. Inverse warp both atlases to native T1

        clear matlabbatch;

        matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = ...
            {yFile};

        matlabbatch{1}.spm.util.defs.comp{1}.inv.space = ...
            {t1};

        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = ...
            {lpba;aal};

        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = ...
            {subDir};

        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = ...
            'native_';

        spm_jobman('run',matlabbatch);

        nativeLPBA = fullfile(subDir, ...
            'native_lpba40.spm5.avg152T1.gm.label.nii');

        nativeAAL = fullfile(subDir, ...
            'native_AAL3v1_1mm.nii');

        assert(isfile(nativeLPBA), ...
            'LPBA40逆变形失败');

        assert(isfile(nativeAAL), ...
            'AAL3逆变形失败');

        %% C. Read CBF

        Vc = spm_vol(cbfCopy);
        C = spm_read_vols(Vc);

        %% D. Resample atlases and GM to CBF grid

        inputs = {nativeLPBA,nativeAAL,gmFile};

        B = cell(3,1);

        [I,J] = ndgrid( ...
            1:Vc.dim(1),1:Vc.dim(2));

        for a = 1:3

            Vs = spm_vol(inputs{a});

            assert(numel(Vs)==1, ...
                '标签或GM不是3D图像');

            M = Vs.mat\Vc.mat;

            B{a} = zeros(Vc.dim);

            interp = 0;
            if a==3
                interp = 1;
            end

            for z = 1:Vc.dim(3)

                X = M(1,1)*I+M(1,2)*J+ ...
                    M(1,3)*z+M(1,4);

                Y = M(2,1)*I+M(2,2)*J+ ...
                    M(2,3)*z+M(2,4);

                Z = M(3,1)*I+M(3,2)*J+ ...
                    M(3,3)*z+M(3,4);

                B{a}(:,:,z) = ...
                    spm_sample_vol(Vs,X,Y,Z,interp);

            end

            B{a}(~isfinite(B{a})) = 0;

            if a<3

                B{a} = round(B{a});
                B{a}(B{a}<0) = 0;

            else

                B{a} = max(0,min(1,B{a}));

            end

        end

        L = B{1};
        A = B{2};
        G = B{3};

        %% E. Save native atlas/GM on CBF grid

        names = { ...
            'LPBA40_on_CBF.nii', ...
            'AAL3_on_CBF.nii', ...
            'GM_on_CBF.nii'};

        for a = 1:3

            Vo = Vc;
            Vo.fname = fullfile(subDir,names{a});
            Vo.pinfo = [1;0;0];

            if a<3
                Vo.dt = [spm_type('int16') 0];
            else
                Vo.dt = [spm_type('float32') 0];
            end

            spm_write_vol(Vo,B{a});

        end

        %% F. LPBA40 55 ROI extraction

        tempL = nan(1,55);
        tempGM_L = zeros(1,55);
        tempFinite_L = zeros(1,55);

        for r = 1:55

            gmMask = (L==LPBA_ID(r)) & ...
                     isfinite(G) & (G>=0.30);

            tempGM_L(r) = nnz(gmMask);

            valid = gmMask & isfinite(C);
            tempFinite_L(r) = nnz(valid);

            if any(valid(:))
                tempL(r) = mean(C(valid));
            end

        end

        %% G. AAL3 150 ROI extraction

        tempA = nan(1,150);
        tempGM_A = zeros(1,150);
        tempFinite_A = zeros(1,150);

        for r = 1:150

            gmMask = (A==AAL_ID(r)) & ...
                     isfinite(G) & (G>=0.30);

            tempGM_A(r) = nnz(gmMask);

            valid = gmMask & isfinite(C);
            tempFinite_A(r) = nnz(valid);

            if any(valid(:))
                tempA(r) = mean(C(valid));
            end

        end

        %% H. Store subject results

        CBF_LPBA(s,:) = tempL;
        CBF_AAL(s,:) = tempA;

        GM_N_LPBA(s,:) = tempGM_L;
        GM_N_AAL(s,:) = tempGM_A;

        Finite_N_LPBA(s,:) = tempFinite_L;
        Finite_N_AAL(s,:) = tempFinite_A;

        success(s) = true;

        fprintf('LPBA40：%d/55 ROI有值\n', ...
            nnz(isfinite(tempL)));

        fprintf('AAL3：%d/150 ROI有值\n', ...
            nnz(isfinite(tempA)));

        if NeedsReview(s)
            fprintf('注意：配准参数较大，建议人工复核\n');
        end

    catch ME

        errorMsg(s) = string(ME.message);

        fprintf('处理失败：%s\n',ME.message);

    end

    %% I. Save checkpoint after every scan

    save(fullfile(outDir,'Progress_Checkpoint.mat'), ...
        'scan','group','pid', ...
        'success','NeedsReview','errorMsg', ...
        'CBF_LPBA','CBF_AAL', ...
        'GM_N_LPBA','GM_N_AAL', ...
        'Finite_N_LPBA','Finite_N_AAL', ...
        'Translation_mm','Rotation_deg', ...
        'LPBA_ID','AAL_ID');

end

%% 11. Build two Excel tables

Base = table( ...
    cellstr(group),cellstr(pid),cellstr(scan), ...
    'VariableNames', ...
    {'Group','Patient_ID','Scan_ID'});

TL = Base;
TA = Base;

for r = 1:55

    TL.(sprintf('ROI_%d',LPBA_ID(r))) = ...
        CBF_LPBA(:,r);

end

for r = 1:150

    TA.(sprintf('ROI_%d',AAL_ID(r))) = ...
        CBF_AAL(:,r);

end

%% 12. Export Excel

lpbaExcel = fullfile(outDir, ...
    'LPBA40_Native_CBF_40Pairs_42H.xlsx');

aalExcel = fullfile(outDir, ...
    'AAL3_Native_CBF_40Pairs_42H.xlsx');

writetable(TL,lpbaExcel,'Sheet','All_CBF');
writetable(TA,aalExcel,'Sheet','All_CBF');

%% 13. Processing log

LPBA_ValidROI = sum(isfinite(CBF_LPBA),2);
AAL_ValidROI = sum(isfinite(CBF_AAL),2);

Log = table( ...
    scan,group,success,NeedsReview, ...
    LPBA_ValidROI,AAL_ValidROI, ...
    errorMsg,Translation_mm,Rotation_deg);

writetable(Log,fullfile(outDir, ...
    'Processing_Log.csv'));

save(fullfile(outDir,'Native_CBF_Results.mat'), ...
    'TL','TA','Log', ...
    'CBF_LPBA','CBF_AAL', ...
    'GM_N_LPBA','GM_N_AAL', ...
    'Finite_N_LPBA','Finite_N_AAL', ...
    'LPBA_ID','AAL_ID');

%% 14. Final summary

fprintf('\n====================================\n');
fprintf('NATIVE SPACE BATCH SUMMARY\n');
fprintf('====================================\n');

fprintf('成功扫描：%d / 122\n',nnz(success));
fprintf('失败扫描：%d\n',nnz(~success));

fprintf('需要检查配准参数：%d\n',nnz(NeedsReview));

fprintf('LPBA40有效CBF：%d / %d\n', ...
    nnz(isfinite(CBF_LPBA)),122*55);

fprintf('AAL3有效CBF：%d / %d\n', ...
    nnz(isfinite(CBF_AAL)),122*150);

if any(~success)

    fprintf('\n失败扫描：\n');
    disp(scan(~success));

end

if any(NeedsReview)

    fprintf('\n需要检查配准的扫描：\n');
    disp(scan(NeedsReview));

end

fprintf('\nLPBA40 Excel：\n%s\n',lpbaExcel);
fprintf('\nAAL3 Excel：\n%s\n',aalExcel);

fprintf('\n全部结果目录：\n%s\n',outDir);

fprintf('\n原始文件未修改。\n');
fprintf('BATCH EXTRACTION COMPLETED.\n');
