/* ==========================================================================
   SAS Analysis Pipeline:
   1. Maternal glucose exposures (FPG & 1-hr PG)
   2. Covariates & gestational weight gain (GWG) classification (WS/T 801-2022)
   3. Anthropometric outcomes (>= 3 visits; WHO z-scores & internal WC z-scores)
   4. Group-based trajectory modeling (PROC TRAJ) grid search & selection
   5. Merging datasets & exporting CSV files
   ========================================================================== */
%let raw_lib_dir  = D:\database\SMBC;
%let data_dir     = D:\database\SMBC\anthrowhodata;
%let out_dir      = D:\database\SMBC\Myresult11r1;
%let traj_inc     = D:\database\SMBC\BMC public health\code data\trajectory.sas;

libname smbc     "&raw_lib_dir";
libname mdt      "&raw_lib_dir\mydata";
libname workdata "&data_dir";

/* ---------- 1. Maternal glucose exposure pathway ---------- */
data mdt.exposure;
    set smbc.total_cmr11;
    if v1 =< 1 then v1 = v1 * 8;
    fpg = v1;

    if v2 ^= . then h1pg2 = v2;
    else if v6 ^= . then h1pg2 = v6;
    else h1pg2 = .;

    if o1 ^= . and c1 ^= . then hlpgw1 = int((o1 - c1) / 7);
    if o2 ^= . and c1 ^= . then hlpgw2 = int((o2 - c1) / 7);
    else if o2 = . and o6 ^= . and c1 ^= . then hlpgw2 = int((o6 - c1) / 7);

    if fpg = . or h1pg2 = . then delete;
    if (12 <= hlpgw1 <= 20) and (20 <= hlpgw2 <= 28) then output;
    keep code hlpgw1 fpg hlpgw2 h1pg2;
run;

/* ---------- 2. Covariates and GWG categorization (WS/T 801-2022) ---------- */
data covars_tmp;
    set smbc.total_cmr11(rename=(
        a4 = mheight a5 = mweight m1 = delierydate
        c1 = lmpdate m4_corr = birthweight A6 = mweight_FV
    ));
    mheightm = mheight / 100;
    mbmi = mweight / (mheightm * mheightm);

    /* Exclusive breastfeeding duration (months) */
    ebf_m = coalesce(sb2, na31d);
    if missing(ebf_m) then ebf_cat = .;
    else if ebf_m < 4 then ebf_cat = 1;
    else if ebf_m < 6 then ebf_cat = 2;
    else ebf_cat = 3;

    /* Continued breastfeeding indicators */
    if sb3 = 1 then bf6 = 1; else if sb3 = 2 then bf6 = 0; else bf6 = .;
    if tb1 = 1 then bf12 = 1; else if tb1 = 2 then bf12 = 0; else bf12 = .;
    if vb1a = 1 then bfv = 1; else if vb1a = 2 then bfv = 0; else bfv = .;

    if missing(bf12) and bfv = 1 then bf12 = 1;
    if missing(bf6) and bf12 = 1 then bf6 = 1;
    if missing(bf12) and bf6 = 0 then bf12 = 0;
    if missing(bfv) and bf12 = 0 then bfv = 0;

    if missing(ebf_cat) then brefeeding = .;
    else if ebf_cat = 3 then brefeeding = 1;
    else brefeeding = 0;

    bweight = birthweight / 1000;
    if . < mbmi < 18.5 then mbmigr = 1;
    else if 18.5 =< mbmi < 24 then mbmigr = 3;
    else if mbmi >= 24 then mbmigr = 4;

    mage = YRDIF(a1, date1, 'AGE');

    if . < m4 < 2500 then lbw = 1;
    else if m4 >= 2500 and m4 < 4000 then lbw = 2;
    else if m4 >= 4000 then lbw = 3;
    if m4 = . then lbw = .;

    /* Total gestational weight gain */
    maxweight = max(of q1-q23);
    if maxweight = . then maxweight = f1a;
    getw = maxweight - mweight;
    if 0 < getw < 10.0 then getw_gr3 = 1;
    else if 10.0 <= getw < 20.0 then getw_gr3 = 2;
    else if getw >= 20.0 then getw_gr3 = 3;

    keep code sex mage parity gweek bweight lbw brefeeding
         fedu fsmoke fwine medu smoke1 income mbmi mbmigr
         mweight maxweight getw getw_gr3 lmpdate delierydate delivery;
run;

data mdt.covars;
    set covars_tmp;
    if not missing(mbmi) and not missing(getw) then do;
        /* Underweight (BMI < 18.5 kg/m2): recommended 11.0-16.0 kg */
        if mbmi < 18.5 then do;
            if getw < 11.0 then gwg_cat = 1;
            else if getw <= 16.0 then gwg_cat = 2;
            else gwg_cat = 3;
        end;
        /* Normal weight (18.5 <= BMI < 24.0 kg/m2): recommended 8.0-14.0 kg */
        else if mbmi < 24.0 then do;
            if getw < 8.0 then gwg_cat = 1;
            else if getw <= 14.0 then gwg_cat = 2;
            else gwg_cat = 3;
        end;
        /* Overweight (24.0 <= BMI < 28.0 kg/m2): recommended 7.0-11.0 kg */
        else if mbmi < 28.0 then do;
            if getw < 7.0 then gwg_cat = 1;
            else if getw <= 11.0 then gwg_cat = 2;
            else gwg_cat = 3;
        end;
        /* Obesity (BMI >= 28.0 kg/m2): recommended 5.0-9.0 kg */
        else do;
            if getw < 5.0 then gwg_cat = 1;
            else if getw <= 9.0 then gwg_cat = 2;
            else gwg_cat = 3;
        end;
    end;
    label gwg_cat = '1=Insufficient, 2=Optimal, 3=Excessive';
run;
proc sort data=mdt.covars; by code; run;

/* ---------- 3. Offspring anthropometric outcome pathway (>= 3 visits) ---------- */
data outcome;
    set smbc.total_cmr11(rename=(
        xf2 = height6m V2XB = height2yr tz6 = height1yr
        ma18b = height4yr ne1a = height6yr M4_corr = wtbirthg
        V2XA = wt2yr TZ5 = wt1yr NE1B = wt6yr ED1B = wt8yr
        weight11 = wt11yr height11 = height11yr
        K3 = wcbirth MT42B = wc4yr NE2B = wc6yr tf5 = wc1yr xf5 = wc6m
    ));
    height8yr = (ed1a1 + ed1a2) / 2;
    wc8yr     = (ed2b1 + ed2b2) / 2;
    height0   = 50;

    if ma18a ^= . and mt41b ^= . then wt4yr = mt41b;
    else if ma18a ^= . and mt41b = . then wt4yr = ma18a;
    else wt4yr = .;

    if xf1 ^= . then wt6m = xf1;
    else wt6m = tz2;

    if wcbirth < 10 then wcbirth = wcbirth * 10;

    wtfreq = n(wt6m, wt1yr, wt2yr, wt4yr, wt6yr, wt8yr, wt11yr);
    htfreq = n(height6m, height1yr, height2yr, height4yr, height6yr, height8yr, height11yr);
    wcfreq = n(wc6m, wc1yr, wc4yr, wc6yr, wc8yr, wc11);
    whfreq = max(wtfreq, htfreq);
    wtbirth = wtbirthg / 1000;

    keep code sex
         wtbirth wt6m wt1yr wt2yr wt4yr wt6yr wt8yr wt11yr
         height0 height6m height1yr height2yr height4yr height6yr height8yr height11yr
         wcbirth wc6m wc1yr wc4yr wc6yr wc8yr wc11 wtfreq htfreq wcfreq whfreq;
run;

data workdata.athropdata;
    set outcome(where=(whfreq >= 3)
                rename=(wt11yr=wt10yr height11yr=height10yr wc11=wc10yr));
    keep code sex
         wtbirth wt6m wt1yr wt2yr wt4yr wt6yr wt8yr wt10yr
         height0 height6m height1yr height2yr height4yr height6yr height8yr height10yr
         wcbirth wc6m wc1yr wc4yr wc6yr wc8yr wc10yr;
run;

/* ---------- 4. WHO external z-scores for weight, height/length, and BMI ---------- */
%include "&data_dir\igrowup-sas\igrowup_standard.sas";

%macro calc_who_0_4(ds_in, age_var, age_val, age_unit, wt_var, ht_var, label_txt);
    data workdata.&ds_in;
        set workdata.athropdata;
        &age_var = &age_val;
        weight = &wt_var;
        height = &ht_var;
        keep code sex &age_var weight height;
    run;
    %igrowup_standard(
        LABEL=&label_txt, REF_LIB=&data_dir\igrowup-sas,
        DATA_LIB=&data_dir, DATA_IN=&ds_in, SEX=sex,
        AGE=&age_var, AGE_UNIT=&age_unit, WEIGHT=weight, LENHEI=height
    );
%mend;

%calc_who_0_4(y_birth, age_days, 0,  DAYS,   wtbirth, height0,   "Birth Z-score");
%calc_who_0_4(y6m,     age_m,    6,  MONTHS, wt6m,    height6m,  "6m");
%calc_who_0_4(y1yr,    age_m,    12, MONTHS, wt1yr,   height1yr, "1yr");
%calc_who_0_4(y2yr,    age_m,    24, MONTHS, wt2yr,   height2yr, "2yr");
%calc_who_0_4(y4yr,    age_m,    48, MONTHS, wt4yr,   height4yr, "4yr");

%include "&data_dir\who2007_sas\who2007.sas";

%macro calc_who_5_10(ds_in, age_val, wt_var, ht_var, label_txt);
    data workdata.&ds_in;
        set workdata.athropdata;
        age_yr = &age_val;
        weight = &wt_var;
        height = &ht_var;
        keep code sex age_yr weight height;
    run;
    %WHO2007(
        LABEL=&label_txt, REF_LIB=&data_dir\who2007_sas,
        DATA_LIB=&data_dir, DATA_IN=&ds_in, SEX=sex,
        AGE=age_yr, AGE_UNIT=YEARS, WEIGHT=weight, HEIGHT=height
    );
%mend;

%calc_who_5_10(y6,  6,  wt6yr,  height6yr,  "6y");
%calc_who_5_10(y8,  8,  wt8yr,  height8yr,  "8y");
%calc_who_5_10(y10, 10, wt10yr, height10yr, "10y");

data z_all_wh;
    merge workdata.athropdata(in=a)
          workdata.y_birth_z_st(keep=code _zwei _zlen _zbmi rename=(_zwei=w1 _zlen=h1 _zbmi=b1))
          workdata.y6m_z_st(keep=code _zwei _zlen _zbmi rename=(_zwei=w2 _zlen=h2 _zbmi=b2))
          workdata.y1yr_z_st(keep=code _zwei _zlen _zbmi rename=(_zwei=w3 _zlen=h3 _zbmi=b3))
          workdata.y2yr_z_st(keep=code _zwei _zlen _zbmi rename=(_zwei=w4 _zlen=h4 _zbmi=b4))
          workdata.y4yr_z_st(keep=code _zwei _zlen _zbmi rename=(_zwei=w5 _zlen=h5 _zbmi=b5))
          workdata.y6_z(keep=code _zwfa _zhfa _zbfa rename=(_zwfa=w6 _zhfa=h6 _zbfa=b6))
          workdata.y8_z(keep=code _zwfa _zhfa _zbfa rename=(_zwfa=w7 _zhfa=h7 _zbfa=b7))
          workdata.y10_z(keep=code _zwfa _zhfa _zbfa rename=(_zwfa=w8 _zhfa=h8 _zbfa=b8));
    by code;
    if a;

    array w_src{8} w1-w8; array h_src{8} h1-h8; array b_src{8} b1-b8;
    array weightz{8} weightz1-weightz8;
    array heightz{8} heightz1-heightz8;
    array bmiz{8} bmiz1-bmiz8;
    array visit{8} visit1-visit8 (0 0.5 1 2 4 6 8 10);

    do i = 1 to 8;
        weightz{i} = ifn(-6 <= w_src{i} <= 5, w_src{i}, .);
        heightz{i} = ifn(-5 <= h_src{i} <= 5, h_src{i}, .);
        bmiz{i}    = ifn(-5 <= b_src{i} <= 5, b_src{i}, .);
    end;
    keep code sex weightz1-weightz8 heightz1-heightz8 bmiz1-bmiz8 visit1-visit8;
run;

/* ---------- 5. Internal sex- and age-standardized waist circumference z-scores ---------- */
data wc_long;
    set workdata.athropdata(rename=(code=id));
    array wc_vals{7} wcbirth wc6m wc1yr wc4yr wc6yr wc8yr wc10yr;
    array time_vals{7} (1 2 3 4 5 6 7);
    array age_vals{7}  (0 0.5 1 4 6 8 10);
    do i = 1 to 7;
        wc = wc_vals{i};
        time_point = time_vals{i};
        visit = age_vals{i};
        if not missing(wc) then output;
    end;
    keep id sex wc time_point visit;
run;

proc sort data=wc_long; by sex visit; run;
proc means data=wc_long noprint;
    by sex visit;
    var wc;
    output out=wc_stats mean=mean_wc std=std_wc;
run;

data wc_long_z;
    merge wc_long wc_stats;
    by sex visit;
    wcz = ifn(std_wc > 0, (wc - mean_wc) / std_wc, 0);
    keep id sex visit time_point wcz;
run;

proc sort data=wc_long_z; by id time_point; run;
data wc_wide_clean;
    set wc_long_z(rename=(id=code));
    by code;
    retain wcz1-wcz7 time1-time7;
    array z_arr{7} wcz1-wcz7;
    array t_arr{7} time1-time7;
    if first.code then call missing(of wcz1-wcz7, of time1-time7);
    z_arr{time_point} = min(max(wcz, -6), 6);
    t_arr{time_point} = visit;
    if last.code then output;
    keep code wcz1-wcz7 time1-time7;
run;

proc sort data=z_all_wh; by code; run;
proc sort data=wc_wide_clean; by code; run;
data workdata.athropdata_8visits_clean;
    merge z_all_wh(in=a) wc_wide_clean;
    by code;
    if a;
run;

/* ---------- 6. Trajectory grid search, model selection, and plotting ---------- */
%macro run_one_outcome(
    code =, cls_prefix =, label_y =, file_prefix =,
    yvars =, tvars =, ymin = -5, ymax = 5,
    min_pct_cut = 15, min_pb_cut = 0.70, outdir = &data_dir
);
    %local ng ord k order_list plist ntot out_ds stat_ds plot_ds
           grp_sum fit_sum fit_fin cand_ds best_ds n_pass best_ng best_ord best_ord_list;

    %let grp_sum = work.&code._grp_summary;
    %let fit_sum = work.&code._fit_summary;
    %let fit_fin = work.&code._fit_final;
    %let cand_ds = work.&code._model_candidates;
    %let best_ds = work.&code._best_model;

    proc datasets lib=work nolist nowarn;
        delete &code._grp_summary &code._fit_summary &code._fit_final
               &code._model_candidates &code._best_model _cls_&cls_prefix;
    quit;

    %do ng = 2 %to 4;
        %do ord = 1 %to 3;
            %let order_list =;
            %do k = 1 %to &ng;
                %let order_list = &order_list &ord;
            %end;

            %let out_ds  = workdata.&code._traj_out_g&ng._o&ord;
            %let stat_ds = workdata.&code._traj_stat_g&ng._o&ord;
            %let plot_ds = workdata.&code._traj_plot_g&ng._o&ord;

            filename trajprt "%sysfunc(pathname(work))/print_&code._g&ng._o&ord..txt";
            proc printto print=trajprt new; run;
            proc traj data=workdata.athropdata_8visits_clean
                      out=&out_ds outstat=&stat_ds outplot=&plot_ds;
                id code; var &yvars; indep &tvars; model cnorm;
                min &ymin; max &ymax; ngroups &ng; order &order_list;
            run;
            proc printto print=print; run;

            data &stat_ds;
                set &stat_ds;
                model_ngroups = &ng;
                model_order   = &ord;
            run;

            /* Parse fit statistics from PROC TRAJ output */
            data work.&code._fit_g&ng._o&ord(keep=ngroups order_used bic_n bic_nstar aic loglik);
                length line $300;
                retain re;
                if _n_ = 1 then
                    re = prxparse('/BIC=\s*(-?\d+\.?\d*)\s*\(N=\s*\d+\)\s*BIC=\s*(-?\d+\.?\d*)\s*\(N=\s*\d+\)\s*AIC=\s*(-?\d+\.?\d*)\s*ll=\s*(-?\d+\.?\d*)/');
                infile trajprt truncover;
                input line $300.;
                ngroups = &ng; order_used = &ord;
                if prxmatch(re, line) then do;
                    bic_n     = input(prxposn(re, 1, line), best32.);
                    bic_nstar = input(prxposn(re, 2, line), best32.);
                    aic       = input(prxposn(re, 3, line), best32.);
                    loglik    = input(prxposn(re, 4, line), best32.);
                    output;
                end;
            run;

            %let plist =;
            %do k = 1 %to &ng; %let plist = &plist grp&k.prb; %end;

            data work._pp_&code;
                set &out_ds;
                array p{&ng} &plist;
                if not missing(group) then pp = p{group};
            run;

            proc sql noprint;
                select count(*) into :ntot trimmed from work._pp_&code where not missing(group);
                create table work._grp_one_&code as
                select &ng as ngroups, &ord as order_used, group,
                       count(*) as n, count(*)/&ntot*100 as pct, mean(pp) as meanprob
                from work._pp_&code where not missing(group) group by group;
            quit;

            proc append base=&grp_sum data=work._grp_one_&code force; run;
        %end;
    %end;

    data &fit_sum;
        set work.&code._fit_g2_o1-work.&code._fit_g4_o3;
    run;
    proc sort data=&fit_sum; by ngroups order_used; run;

    proc sql;
        create table &fit_fin as
        select ngroups, order_used, bic_n, bic_nstar, aic, loglik
        from &fit_sum group by ngroups, order_used having loglik = max(loglik);

        create table &cand_ds as
        select a.ngroups, a.order_used, a.bic_nstar, a.aic, a.loglik,
               b.min_pct, b.min_meanprob,
               case when b.min_meanprob >= &min_pb_cut and b.min_pct >= &min_pct_cut
                    then "Pass" else "Fail" end as pass_23
        from &fit_fin as a
        inner join (
            select ngroups, order_used, min(pct) as min_pct, min(meanprob) as min_meanprob
            from &grp_sum group by ngroups, order_used
        ) as b
        on a.ngroups = b.ngroups and a.order_used = b.order_used
        order by bic_nstar desc, aic desc, loglik desc;
    quit;

    %let n_pass = 0;
    proc sql noprint;
        select count(*) into :n_pass trimmed from &cand_ds where pass_23 = "Pass";
    quit;

    %if &n_pass > 0 %then %do;
        data &best_ds;
            set &cand_ds(where=(pass_23 = "Pass") obs=1);
        run;
        proc sql noprint;
            select ngroups, order_used into :best_ng trimmed, :best_ord trimmed from &best_ds;
            select count(*) into :ntot trimmed from workdata.&code._traj_out_g&best_ng._o&best_ord where not missing(group);
            create table work._lab_&code as
            select group, count(*)/&ntot*100 as pct
            from workdata.&code._traj_out_g&best_ng._o&best_ord where not missing(group) group by group;
        quit;

        data _null_;
            set work._lab_&code;
            length grp_name $40;
            if &best_ng = 2 then do;
                if group = 1 then grp_name = "Stable group";
                else if group = 2 then grp_name = "Increased group";
                else grp_name = cats("Group ", group);
            end;
            else do;
                if group = 1 then grp_name = "Stable group";
                else if group = 2 then grp_name = "Moderate increased group";
                else if group = 3 then grp_name = "Sharp increase group";
                else grp_name = cats("Group ", group);
            end;
            call symputx(cats("lab_&code._", group), cats(grp_name, " (", put(pct, 4.1), "%)"));
        run;

        data work.plot_&code._g&best_ng._o&best_ord._long;
            set workdata.&code._traj_plot_g&best_ng._o&best_ord;
            length grp $55;
            array avgs{&best_ng}  avg1-avg&best_ng;
            array preds{&best_ng} pred1-pred&best_ng;
            array lows{&best_ng}  l95m1-l95m&best_ng;
            array highs{&best_ng} u95m1-u95m&best_ng;
            do k = 1 to &best_ng;
                grp = symget(cats("lab_&code._", k));
                avg = avgs{k}; pred = preds{k}; low = lows{k}; high = highs{k};
                output;
            end;
            keep t grp avg pred low high;
        run;

        ods html close;
        ods listing gpath="&outdir" image_dpi=300;
        ods graphics on / reset=all imagename="&file_prefix._g&best_ng._o&best_ord"
                         imagefmt=tiff width=9in height=6in;
        proc sgplot data=work.plot_&code._g&best_ng._o&best_ord._long noborder;
            band x=t lower=low upper=high / group=grp transparency=0.7;
            series x=t y=pred / group=grp lineattrs=(thickness=2) name="s";
            scatter x=t y=avg / group=grp markerattrs=(symbol=circlefilled size=6);
            xaxis label="Follow-up time (years)" labelattrs=(size=12) valueattrs=(size=10);
            yaxis label="&label_y" labelattrs=(size=12) valueattrs=(size=10);
            keylegend "s" / title="" valueattrs=(size=8.5) linelength=15
                            across=&best_ng noborder location=outside position=bottom;
        run;
        ods graphics off; ods listing close; ods html;

        %let best_ord_list =;
        %do k = 1 %to &best_ng; %let best_ord_list = &best_ord_list &best_ord; %end;

        proc traj data=workdata.athropdata_8visits_clean
                  out=work._best_o_&code outstat=work._best_s_&code outplot=work._best_p_&code;
            id code; var &yvars; indep &tvars; model cnorm;
            min &ymin; max &ymax; ngroups &best_ng; order &best_ord_list;
        run;

        data work._cls_&cls_prefix;
            set work._best_o_&code(rename=(group=&cls_prefix._class));
            keep code &cls_prefix._class;
        run;
        proc sort data=work._cls_&cls_prefix; by code; run;
    %end;
%mend;

%run_one_outcome(code=wt,  cls_prefix=weightz, label_y=Weight Z-score, file_prefix=Weight_traj, yvars=weightz1-weightz8, tvars=visit1-visit8, ymin=-5,  ymax=5);
%run_one_outcome(code=ht,  cls_prefix=heightz, label_y=Height Z-score, file_prefix=Height_traj, yvars=heightz1-heightz8, tvars=visit1-visit8, ymin=-5,  ymax=5);
%run_one_outcome(code=bmi, cls_prefix=bmiz,    label_y=BMI Z-score,    file_prefix=BMI_traj,    yvars=bmiz1-bmiz8,       tvars=visit1-visit8, ymin=-15, ymax=15);
%run_one_outcome(code=wcz, cls_prefix=wcz,     label_y=WC Z-score,     file_prefix=WC_traj,     yvars=wcz1-wcz7,         tvars=time1-time7,   ymin=-6,  ymax=6);

data workdata.zscore_class_bymode;
    merge work._cls_bmiz(in=in_bmi) work._cls_weightz work._cls_heightz work._cls_wcz;
    by code;
    if in_bmi;
run;

/* ---------- 7. Linkage across pathways and export final datasets ---------- */
%include "&traj_inc";

data traj_all;
    update workdata.trajectory_data workdata.zscore_class_bymode;
    by code;
run;

data traj_renamed;
    set traj_all(rename=(
        bmiz_class=bmiz_traj_model
        weightz_class=wtz_traj_model
        heightz_class=htz_traj_model
        wcz_class=wcz_traj_model
    ));
run;
proc sort data=traj_renamed; by code; run;
proc sort data=mdt.exposure out=exp_sorted; by code; run;

data merged_555;
    merge exp_sorted(in=a) traj_renamed(in=b) mdt.covars(in=c) workdata.owo_who_20261113(in=d);
    by code;
    if a and b and c and d;
run;

proc rank data=merged_555 out=ranked_3 groups=3;
    var fpg h1pg2;
    ranks fpg_gr3 pg1h_gr3;
run;

proc rank data=ranked_3 out=ranked_5 groups=5;
    var fpg h1pg2;
    ranks fpg_gr5 pg1h_gr5;
run;

data workdata.dataset261115;
    set ranked_5;
    array grp_vars{4} fpg_gr3 pg1h_gr3 fpg_gr5 pg1h_gr5;
    do i = 1 to 4;
        if not missing(grp_vars{i}) then grp_vars{i} = grp_vars{i} + 1;
    end;
    drop i;
run;

proc sql;
    create table allinex as
    select c.*, case when e.code is not null then 1 else 0 end as in_exposure
    from mdt.covars as c
    left join workdata.dataset261115 as e
    on c.code = e.code;
quit;

options dlcreatedir;
libname exp_out "&out_dir";
proc export data=allinex outfile="&out_dir\inandex1115.csv" dbms=csv replace; run;
proc export data=workdata.dataset261115 outfile="&out_dir\tertile16.csv" dbms=csv replace; run;
