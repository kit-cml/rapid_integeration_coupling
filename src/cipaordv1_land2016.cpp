/*
   There are a total of 200 entries in the algebraic variable array.
   There are a total of 49 entries in each of the rate and state variable arrays.
   There are a total of 206 entries in the constant variable array.
 */

#include "cipaordv1_land2016.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

// #include <Eigen/Dense>
// #include <Eigen/Eigenvalues>

// using Eigen::MatrixXd;
// using Eigen::MatrixXcd;
// using Eigen::VectorXcd;
// using Eigen::EigenSolver;

/*
 * TIME is time in component environment (millisecond).
 * CONSTANTS[cipaordv1_celltype] is celltype in component environment (dimensionless).
 * CONSTANTS[cipaordv1_nao] is nao in component extracellular (millimolar).
 * CONSTANTS[cipaordv1_cao] is cao in component extracellular (millimolar).
 * CONSTANTS[cipaordv1_ko] is ko in component extracellular (millimolar).
 * CONSTANTS[cipaordv1_R] is R in component physical_constants (joule_per_kilomole_kelvin).
 * CONSTANTS[cipaordv1_T] is T in component physical_constants (kelvin).
 * CONSTANTS[cipaordv1_F] is F in component physical_constants (coulomb_per_mole).
 * CONSTANTS[cipaordv1_zna] is zna in component physical_constants (dimensionless).
 * CONSTANTS[cipaordv1_zca] is zca in component physical_constants (dimensionless).
 * CONSTANTS[cipaordv1_zk] is zk in component physical_constants (dimensionless).
 * CONSTANTS[cipaordv1_L] is L in component cell_geometry (centimeter).
 * CONSTANTS[cipaordv1_rad] is rad in component cell_geometry (centimeter).
 * CONSTANTS[cipaordv1_vcell] is vcell in component cell_geometry (microliter).
 * CONSTANTS[cipaordv1_Ageo] is Ageo in component cell_geometry (centimeter_squared).
 * CONSTANTS[cipaordv1_Acap] is Acap in component cell_geometry (centimeter_squared).
 * CONSTANTS[cipaordv1_vmyo] is vmyo in component cell_geometry (microliter).
 * CONSTANTS[cipaordv1_vnsr] is vnsr in component cell_geometry (microliter).
 * CONSTANTS[cipaordv1_vjsr] is vjsr in component cell_geometry (microliter).
 * CONSTANTS[cipaordv1_vss] is vss in component cell_geometry (microliter).
 * STATES[cipaordv1_V] is v in component membrane (millivolt).
 * ALGEBRAIC[cipaordv1_vfrt] is vfrt in component membrane (dimensionless).
 * CONSTANTS[cipaordv1_ffrt] is ffrt in component membrane (coulomb_per_mole_millivolt).
 * CONSTANTS[cipaordv1_frt] is frt in component membrane (per_millivolt).
 * ALGEBRAIC[cipaordv1_INa] is INa in component INa (microA_per_microF).
 * ALGEBRAIC[cipaordv1_INaL] is INaL in component INaL (microA_per_microF).
 * ALGEBRAIC[cipaordv1_Ito] is Ito in component Ito (microA_per_microF).
 * ALGEBRAIC[cipaordv1_ICaL] is ICaL in component ICaL (microA_per_microF).
 * ALGEBRAIC[cipaordv1_ICaNa] is ICaNa in component ICaL (microA_per_microF).
 * ALGEBRAIC[cipaordv1_ICaK] is ICaK in component ICaL (microA_per_microF).
 * ALGEBRAIC[cipaordv1_IKr] is IKr in component IKr (microA_per_microF).
 * ALGEBRAIC[cipaordv1_IKs] is IKs in component IKs (microA_per_microF).
 * ALGEBRAIC[cipaordv1_IK1] is IK1 in component IK1 (microA_per_microF).
 * ALGEBRAIC[cipaordv1_INaCa_i] is INaCa_i in component INaCa_i (microA_per_microF).
 * ALGEBRAIC[cipaordv1_INaCa_ss] is INaCa_ss in component INaCa_i (microA_per_microF).
 * ALGEBRAIC[cipaordv1_INaK] is INaK in component INaK (microA_per_microF).
 * ALGEBRAIC[cipaordv1_INab] is INab in component INab (microA_per_microF).
 * ALGEBRAIC[cipaordv1_IKb] is IKb in component IKb (microA_per_microF).
 * ALGEBRAIC[cipaordv1_IpCa] is IpCa in component IpCa (microA_per_microF).
 * ALGEBRAIC[cipaordv1_ICab] is ICab in component ICab (microA_per_microF).
 * ALGEBRAIC[cipaordv1_Istim] is Istim in component membrane (microA_per_microF).
 * CONSTANTS[cipaordv1_stim_start] is stim_start in component membrane (millisecond).
 * CONSTANTS[cipaordv1_stim_end] is stim_end in component membrane (millisecond).
 * CONSTANTS[cipaordv1_amp] is amp in component membrane (microA_per_microF).
 * CONSTANTS[cipaordv1_BCL] is BCL in component membrane (millisecond).
 * CONSTANTS[cipaordv1_duration] is duration in component membrane (millisecond).
 * CONSTANTS[cipaordv1_KmCaMK] is KmCaMK in component CaMK (millimolar).
 * CONSTANTS[cipaordv1_aCaMK] is aCaMK in component CaMK (per_millimolar_per_millisecond).
 * CONSTANTS[cipaordv1_bCaMK] is bCaMK in component CaMK (per_millisecond).
 * CONSTANTS[cipaordv1_CaMKo] is CaMKo in component CaMK (dimensionless).
 * CONSTANTS[cipaordv1_KmCaM] is KmCaM in component CaMK (millimolar).
 * ALGEBRAIC[cipaordv1_CaMKb] is CaMKb in component CaMK (millimolar).
 * ALGEBRAIC[cipaordv1_CaMKa] is CaMKa in component CaMK (millimolar).
 * STATES[cipaordv1_CaMKt] is CaMKt in component CaMK (millimolar).
 * STATES[cipaordv1_cass] is cass in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_cmdnmax_b] is cmdnmax_b in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_cmdnmax] is cmdnmax in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_kmcmdn] is kmcmdn in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_trpnmax] is trpnmax in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_kmtrpn] is kmtrpn in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_BSRmax] is BSRmax in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_KmBSR] is KmBSR in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_BSLmax] is BSLmax in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_KmBSL] is KmBSL in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_csqnmax] is csqnmax in component intracellular_ions (millimolar).
 * CONSTANTS[cipaordv1_kmcsqn] is kmcsqn in component intracellular_ions (millimolar).
 * STATES[cipaordv1_nai] is nai in component intracellular_ions (millimolar).
 * STATES[cipaordv1_nass] is nass in component intracellular_ions (millimolar).
 * STATES[cipaordv1_ki] is ki in component intracellular_ions (millimolar).
 * STATES[cipaordv1_kss] is kss in component intracellular_ions (millimolar).
 * STATES[cipaordv1_cansr] is cansr in component intracellular_ions (millimolar).
 * STATES[cipaordv1_cajsr] is cajsr in component intracellular_ions (millimolar).
 * STATES[cipaordv1_cai] is cai in component intracellular_ions (millimolar).
 * ALGEBRAIC[cipaordv1_JdiffNa] is JdiffNa in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_Jdiff] is Jdiff in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_Jup] is Jup in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_JdiffK] is JdiffK in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_Jrel] is Jrel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_Jtr] is Jtr in component trans_flux (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_Bcai] is Bcai in component intracellular_ions (dimensionless).
 * ALGEBRAIC[cipaordv1_Bcajsr] is Bcajsr in component intracellular_ions (dimensionless).
 * ALGEBRAIC[cipaordv1_Bcass] is Bcass in component intracellular_ions (dimensionless).
 * CONSTANTS[cipaordv1_cm] is cm in component intracellular_ions (microF_per_centimeter_squared).
 * CONSTANTS[cipaordv1_PKNa] is PKNa in component reversal_potentials (dimensionless).
 * ALGEBRAIC[cipaordv1_ENa] is ENa in component reversal_potentials (millivolt).
 * ALGEBRAIC[cipaordv1_EK] is EK in component reversal_potentials (millivolt).
 * ALGEBRAIC[cipaordv1_EKs] is EKs in component reversal_potentials (millivolt).
 * ALGEBRAIC[cipaordv1_mss] is mss in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_tm] is tm in component INa (millisecond).
 * CONSTANTS[cipaordv1_mssV1] is mssV1 in component INa (millivolt).
 * CONSTANTS[cipaordv1_mssV2] is mssV2 in component INa (millivolt).
 * CONSTANTS[cipaordv1_mtV1] is mtV1 in component INa (millivolt).
 * CONSTANTS[cipaordv1_mtV2] is mtV2 in component INa (millivolt).
 * CONSTANTS[cipaordv1_mtD1] is mtD1 in component INa (dimensionless).
 * CONSTANTS[cipaordv1_mtD2] is mtD2 in component INa (dimensionless).
 * CONSTANTS[cipaordv1_mtV3] is mtV3 in component INa (millivolt).
 * CONSTANTS[cipaordv1_mtV4] is mtV4 in component INa (millivolt).
 * STATES[cipaordv1_m] is m in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_hss] is hss in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_thf] is thf in component INa (millisecond).
 * ALGEBRAIC[cipaordv1_ths] is ths in component INa (millisecond).
 * CONSTANTS[cipaordv1_hssV1] is hssV1 in component INa (millivolt).
 * CONSTANTS[cipaordv1_hssV2] is hssV2 in component INa (millivolt).
 * CONSTANTS[cipaordv1_Ahs] is Ahs in component INa (dimensionless).
 * CONSTANTS[cipaordv1_Ahf] is Ahf in component INa (dimensionless).
 * STATES[cipaordv1_hf] is hf in component INa (dimensionless).
 * STATES[cipaordv1_hs] is hs in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_h] is h in component INa (dimensionless).
 * CONSTANTS[cipaordv1_GNa] is GNa in component INa (milliS_per_microF).
 * CONSTANTS[cipaordv1_shift_INa_inact] is shift_INa_inact in component INa (millivolt).
 * ALGEBRAIC[cipaordv1_jss] is jss in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_tj] is tj in component INa (millisecond).
 * STATES[cipaordv1_j] is j in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_hssp] is hssp in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_thsp] is thsp in component INa (millisecond).
 * STATES[cipaordv1_hsp] is hsp in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_hp] is hp in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_tjp] is tjp in component INa (millisecond).
 * STATES[cipaordv1_jp] is jp in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_fINap] is fINap in component INa (dimensionless).
 * ALGEBRAIC[cipaordv1_mLss] is mLss in component INaL (dimensionless).
 * ALGEBRAIC[cipaordv1_tmL] is tmL in component INaL (millisecond).
 * STATES[cipaordv1_mL] is mL in component INaL (dimensionless).
 * CONSTANTS[cipaordv1_thL] is thL in component INaL (millisecond).
 * ALGEBRAIC[cipaordv1_hLss] is hLss in component INaL (dimensionless).
 * STATES[cipaordv1_hL] is hL in component INaL (dimensionless).
 * ALGEBRAIC[cipaordv1_hLssp] is hLssp in component INaL (dimensionless).
 * CONSTANTS[cipaordv1_thLp] is thLp in component INaL (millisecond).
 * STATES[cipaordv1_hLp] is hLp in component INaL (dimensionless).
 * CONSTANTS[cipaordv1_GNaL_b] is GNaL_b in component INaL (milliS_per_microF).
 * CONSTANTS[cipaordv1_GNaL] is GNaL in component INaL (milliS_per_microF).
 * ALGEBRAIC[cipaordv1_fINaLp] is fINaLp in component INaL (dimensionless).
 * CONSTANTS[cipaordv1_Gto_b] is Gto_b in component Ito (milliS_per_microF).
 * ALGEBRAIC[cipaordv1_ass] is ass in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_ta] is ta in component Ito (millisecond).
 * STATES[cipaordv1_a] is a in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_iss] is iss in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_delta_epi] is delta_epi in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_tiF_b] is tiF_b in component Ito (millisecond).
 * ALGEBRAIC[cipaordv1_tiS_b] is tiS_b in component Ito (millisecond).
 * ALGEBRAIC[cipaordv1_tiF] is tiF in component Ito (millisecond).
 * ALGEBRAIC[cipaordv1_tiS] is tiS in component Ito (millisecond).
 * ALGEBRAIC[cipaordv1_AiF] is AiF in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_AiS] is AiS in component Ito (dimensionless).
 * STATES[cipaordv1_iF] is iF in component Ito (dimensionless).
 * STATES[cipaordv1_iS] is iS in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_i] is i in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_assp] is assp in component Ito (dimensionless).
 * STATES[cipaordv1_ap] is ap in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_dti_develop] is dti_develop in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_dti_recover] is dti_recover in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_tiFp] is tiFp in component Ito (millisecond).
 * ALGEBRAIC[cipaordv1_tiSp] is tiSp in component Ito (millisecond).
 * STATES[cipaordv1_iFp] is iFp in component Ito (dimensionless).
 * STATES[cipaordv1_iSp] is iSp in component Ito (dimensionless).
 * ALGEBRAIC[cipaordv1_ip] is ip in component Ito (dimensionless).
 * CONSTANTS[cipaordv1_Gto] is Gto in component Ito (milliS_per_microF).
 * ALGEBRAIC[cipaordv1_fItop] is fItop in component Ito (dimensionless).
 * CONSTANTS[cipaordv1_Kmn] is Kmn in component ICaL (millimolar).
 * CONSTANTS[cipaordv1_k2n] is k2n in component ICaL (per_millisecond).
 * CONSTANTS[cipaordv1_PCa_b] is PCa_b in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_dss] is dss in component ICaL (dimensionless).
 * STATES[cipaordv1_d] is d in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_fss] is fss in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_Aff] is Aff in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_Afs] is Afs in component ICaL (dimensionless).
 * STATES[cipaordv1_ff] is ff in component ICaL (dimensionless).
 * STATES[cipaordv1_fs] is fs in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_f] is f in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_fcass] is fcass in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_Afcaf] is Afcaf in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_Afcas] is Afcas in component ICaL (dimensionless).
 * STATES[cipaordv1_fcaf] is fcaf in component ICaL (dimensionless).
 * STATES[cipaordv1_fcas] is fcas in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_fca] is fca in component ICaL (dimensionless).
 * STATES[cipaordv1_jca] is jca in component ICaL (dimensionless).
 * STATES[cipaordv1_ffp] is ffp in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_fp] is fp in component ICaL (dimensionless).
 * STATES[cipaordv1_fcafp] is fcafp in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_fcap] is fcap in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_km2n] is km2n in component ICaL (per_millisecond).
 * ALGEBRAIC[cipaordv1_anca] is anca in component ICaL (dimensionless).
 * STATES[cipaordv1_nca] is nca in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_PhiCaL] is PhiCaL in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_PhiCaNa] is PhiCaNa in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_PhiCaK] is PhiCaK in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_PCa] is PCa in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_PCap] is PCap in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_PCaNa] is PCaNa in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_PCaK] is PCaK in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_PCaNap] is PCaNap in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_PCaKp] is PCaKp in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_fICaLp] is fICaLp in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_td] is td in component ICaL (millisecond).
 * ALGEBRAIC[cipaordv1_tff] is tff in component ICaL (millisecond).
 * ALGEBRAIC[cipaordv1_tfs] is tfs in component ICaL (millisecond).
 * ALGEBRAIC[cipaordv1_tfcaf] is tfcaf in component ICaL (millisecond).
 * ALGEBRAIC[cipaordv1_tfcas] is tfcas in component ICaL (millisecond).
 * CONSTANTS[cipaordv1_tjca] is tjca in component ICaL (millisecond).
 * ALGEBRAIC[cipaordv1_tffp] is tffp in component ICaL (millisecond).
 * ALGEBRAIC[cipaordv1_tfcafp] is tfcafp in component ICaL (millisecond).
 * CONSTANTS[cipaordv1_v0_CaL] is v0 in component ICaL (millivolt).
 * ALGEBRAIC[cipaordv1_A_1] is A_1 in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_B_1] is B_1 in component ICaL (per_millivolt).
 * ALGEBRAIC[cipaordv1_U_1] is U_1 in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_A_2] is A_2 in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_B_2] is B_2 in component ICaL (per_millivolt).
 * ALGEBRAIC[cipaordv1_U_2] is U_2 in component ICaL (dimensionless).
 * ALGEBRAIC[cipaordv1_A_3] is A_3 in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_B_3] is B_3 in component ICaL (per_millivolt).
 * ALGEBRAIC[cipaordv1_U_3] is U_3 in component ICaL (dimensionless).
 * CONSTANTS[cipaordv1_GKr_b] is GKr_b in component IKr (milliS_per_microF).
 * STATES[cipaordv1_IC1] is IC1 in component IKr (dimensionless).
 * STATES[cipaordv1_IC2] is IC2 in component IKr (dimensionless).
 * STATES[cipaordv1_C1] is C1 in component IKr (dimensionless).
 * STATES[cipaordv1_C2] is C2 in component IKr (dimensionless).
 * STATES[cipaordv1_O] is O in component IKr (dimensionless).
 * STATES[cipaordv1_IO] is IO in component IKr (dimensionless).
 * STATES[cipaordv1_IObound] is IObound in component IKr (dimensionless).
 * STATES[cipaordv1_Obound] is Obound in component IKr (dimensionless).
 * STATES[cipaordv1_Cbound] is Cbound in component IKr (dimensionless).
 * STATES[cipaordv1_D] is D in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_GKr] is GKr in component IKr (milliS_per_microF).
 * CONSTANTS[cipaordv1_A1] is A1 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B1] is B1 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q1] is q1 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A2] is A2 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B2] is B2 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q2] is q2 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A3] is A3 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B3] is B3 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q3] is q3 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A4] is A4 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B4] is B4 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q4] is q4 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A11] is A11 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B11] is B11 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q11] is q11 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A21] is A21 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B21] is B21 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q21] is q21 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A31] is A31 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B31] is B31 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q31] is q31 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A41] is A41 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B41] is B41 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q41] is q41 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A51] is A51 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B51] is B51 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q51] is q51 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A52] is A52 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B52] is B52 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q52] is q52 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A53] is A53 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B53] is B53 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q53] is q53 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A61] is A61 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B61] is B61 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q61] is q61 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A62] is A62 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B62] is B62 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q62] is q62 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_A63] is A63 in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_B63] is B63 in component IKr (per_millivolt).
 * CONSTANTS[cipaordv1_q63] is q63 in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_Kmax] is Kmax in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_Ku] is Ku in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_n] is n in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_halfmax] is halfmax in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_Kt] is Kt in component IKr (per_millisecond).
 * CONSTANTS[cipaordv1_Vhalf] is Vhalf in component IKr (millivolt).
 * CONSTANTS[cipaordv1_Temp] is Temp in component IKr (dimensionless).
 * CONSTANTS[cipaordv1_GKs_b] is GKs_b in component IKs (milliS_per_microF).
 * CONSTANTS[cipaordv1_GKs] is GKs in component IKs (milliS_per_microF).
 * ALGEBRAIC[cipaordv1_xs1ss] is xs1ss in component IKs (dimensionless).
 * ALGEBRAIC[cipaordv1_xs2ss] is xs2ss in component IKs (dimensionless).
 * ALGEBRAIC[cipaordv1_txs1] is txs1 in component IKs (millisecond).
 * CONSTANTS[cipaordv1_txs1_max] is txs1_max in component IKs (millisecond).
 * STATES[cipaordv1_xs1] is xs1 in component IKs (dimensionless).
 * STATES[cipaordv1_xs2] is xs2 in component IKs (dimensionless).
 * ALGEBRAIC[cipaordv1_KsCa] is KsCa in component IKs (dimensionless).
 * ALGEBRAIC[cipaordv1_txs2] is txs2 in component IKs (millisecond).
 * CONSTANTS[cipaordv1_GK1] is GK1 in component IK1 (milliS_per_microF).
 * CONSTANTS[cipaordv1_GK1_b] is GK1_b in component IK1 (milliS_per_microF).
 * ALGEBRAIC[cipaordv1_xk1ss] is xk1ss in component IK1 (dimensionless).
 * ALGEBRAIC[cipaordv1_txk1] is txk1 in component IK1 (millisecond).
 * STATES[cipaordv1_xk1] is xk1 in component IK1 (dimensionless).
 * ALGEBRAIC[cipaordv1_rk1] is rk1 in component IK1 (millisecond).
 * CONSTANTS[cipaordv1_kna1] is kna1 in component INaCa_i (per_millisecond).
 * CONSTANTS[cipaordv1_kna2] is kna2 in component INaCa_i (per_millisecond).
 * CONSTANTS[cipaordv1_kna3] is kna3 in component INaCa_i (per_millisecond).
 * CONSTANTS[cipaordv1_kasymm] is kasymm in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_wna] is wna in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_wca] is wca in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_wnaca] is wnaca in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_kcaon] is kcaon in component INaCa_i (per_millisecond).
 * CONSTANTS[cipaordv1_kcaoff] is kcaoff in component INaCa_i (per_millisecond).
 * CONSTANTS[cipaordv1_qna] is qna in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_qca] is qca in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_hna] is hna in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_hca] is hca in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_KmCaAct] is KmCaAct in component INaCa_i (millimolar).
 * CONSTANTS[cipaordv1_Gncx_b] is Gncx_b in component INaCa_i (milliS_per_microF).
 * CONSTANTS[cipaordv1_Gncx] is Gncx in component INaCa_i (milliS_per_microF).
 * ALGEBRAIC[cipaordv1_h1_i] is h1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h2_i] is h2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h3_i] is h3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h4_i] is h4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h5_i] is h5_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h6_i] is h6_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h7_i] is h7_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h8_i] is h8_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h9_i] is h9_i in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_h10_i] is h10_i in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_h11_i] is h11_i in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_h12_i] is h12_i in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_k1_i] is k1_i in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_k2_i] is k2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k3p_i] is k3p_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k3pp_i] is k3pp_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k3_i] is k3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k4_i] is k4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k4p_i] is k4p_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k4pp_i] is k4pp_i in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_k5_i] is k5_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k6_i] is k6_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k7_i] is k7_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k8_i] is k8_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_x1_i] is x1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_x2_i] is x2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_x3_i] is x3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_x4_i] is x4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_E1_i] is E1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_E2_i] is E2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_E3_i] is E3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_E4_i] is E4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_allo_i] is allo_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_JncxNa_i] is JncxNa_i in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_JncxCa_i] is JncxCa_i in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_h1_ss] is h1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h2_ss] is h2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h3_ss] is h3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h4_ss] is h4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h5_ss] is h5_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h6_ss] is h6_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h7_ss] is h7_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h8_ss] is h8_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_h9_ss] is h9_ss in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_h10_ss] is h10_ss in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_h11_ss] is h11_ss in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_h12_ss] is h12_ss in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_k1_ss] is k1_ss in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_k2_ss] is k2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k3p_ss] is k3p_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k3pp_ss] is k3pp_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k3_ss] is k3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k4_ss] is k4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k4p_ss] is k4p_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k4pp_ss] is k4pp_ss in component INaCa_i (dimensionless).
 * CONSTANTS[cipaordv1_k5_ss] is k5_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k6_ss] is k6_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k7_ss] is k7_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_k8_ss] is k8_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_x1_ss] is x1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_x2_ss] is x2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_x3_ss] is x3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_x4_ss] is x4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_E1_ss] is E1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_E2_ss] is E2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_E3_ss] is E3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_E4_ss] is E4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_allo_ss] is allo_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[cipaordv1_JncxNa_ss] is JncxNa_ss in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_JncxCa_ss] is JncxCa_ss in component INaCa_i (millimolar_per_millisecond).
 * CONSTANTS[cipaordv1_k1p] is k1p in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_k1m] is k1m in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_k2p] is k2p in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_k2m] is k2m in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_k3p] is k3p in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_k3m] is k3m in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_k4p] is k4p in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_k4m] is k4m in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_Knai0] is Knai0 in component INaK (millimolar).
 * CONSTANTS[cipaordv1_Knao0] is Knao0 in component INaK (millimolar).
 * CONSTANTS[cipaordv1_delta] is delta in component INaK (millivolt).
 * CONSTANTS[cipaordv1_Kki] is Kki in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_Kko] is Kko in component INaK (per_millisecond).
 * CONSTANTS[cipaordv1_MgADP] is MgADP in component INaK (millimolar).
 * CONSTANTS[cipaordv1_MgATP] is MgATP in component INaK (millimolar).
 * CONSTANTS[cipaordv1_Kmgatp] is Kmgatp in component INaK (millimolar).
 * CONSTANTS[cipaordv1_H] is H in component INaK (millimolar).
 * CONSTANTS[cipaordv1_eP] is eP in component INaK (dimensionless).
 * CONSTANTS[cipaordv1_Khp] is Khp in component INaK (millimolar).
 * CONSTANTS[cipaordv1_Knap] is Knap in component INaK (millimolar).
 * CONSTANTS[cipaordv1_Kxkur] is Kxkur in component INaK (millimolar).
 * CONSTANTS[cipaordv1_Pnak_b] is Pnak_b in component INaK (milliS_per_microF).
 * CONSTANTS[cipaordv1_Pnak] is Pnak in component INaK (milliS_per_microF).
 * ALGEBRAIC[cipaordv1_Knai] is Knai in component INaK (millimolar).
 * ALGEBRAIC[cipaordv1_Knao] is Knao in component INaK (millimolar).
 * ALGEBRAIC[cipaordv1_P] is P in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_a1] is a1 in component INaK (dimensionless).
 * CONSTANTS[cipaordv1_b1] is b1 in component INaK (dimensionless).
 * CONSTANTS[cipaordv1_a2] is a2 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_b2] is b2 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_a3] is a3 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_b3] is b3 in component INaK (dimensionless).
 * CONSTANTS[cipaordv1_a4] is a4 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_b4] is b4 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_x1] is x1 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_x2] is x2 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_x3] is x3 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_x4] is x4 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_E1] is E1 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_E2] is E2 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_E3] is E3 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_E4] is E4 in component INaK (dimensionless).
 * ALGEBRAIC[cipaordv1_JnakNa] is JnakNa in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_JnakK] is JnakK in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_xkb] is xkb in component IKb (dimensionless).
 * CONSTANTS[cipaordv1_GKb_b] is GKb_b in component IKb (milliS_per_microF).
 * CONSTANTS[cipaordv1_GKb] is GKb in component IKb (milliS_per_microF).
 * CONSTANTS[cipaordv1_PNab] is PNab in component INab (milliS_per_microF).
 * ALGEBRAIC[cipaordv1_A_Nab] is A in component INab (microA_per_microF).
 * CONSTANTS[cipaordv1_B_Nab] is B in component INab (per_millivolt).
 * CONSTANTS[cipaordv1_v0_Nab] is v0 in component INab (millivolt).
 * ALGEBRAIC[cipaordv1_U] is U in component INab (dimensionless).
 * CONSTANTS[cipaordv1_PCab] is PCab in component ICab (milliS_per_microF).
 * ALGEBRAIC[cipaordv1_A_Cab] is A in component ICab (microA_per_microF).
 * CONSTANTS[cipaordv1_B_Cab] is B in component ICab (per_millivolt).
 * CONSTANTS[cipaordv1_v0_Cab] is v0 in component ICab (millivolt).
 * ALGEBRAIC[cipaordv1_U] is U in component ICab (dimensionless).
 * CONSTANTS[cipaordv1_GpCa] is GpCa in component IpCa (milliS_per_microF).
 * CONSTANTS[cipaordv1_KmCap] is KmCap in component IpCa (millimolar).
 * CONSTANTS[cipaordv1_bt] is bt in component ryr (millisecond).
 * CONSTANTS[cipaordv1_a_rel] is a_rel in component ryr (millisecond).
 * ALGEBRAIC[cipaordv1_Jrel_inf] is Jrel_inf in component ryr (dimensionless).
 * ALGEBRAIC[cipaordv1_tau_rel] is tau_rel in component ryr (millisecond).
 * ALGEBRAIC[cipaordv1_Jrel_infp] is Jrel_infp in component ryr (dimensionless).
 * ALGEBRAIC[cipaordv1_Jrel_temp] is Jrel_temp in component ryr (dimensionless).
 * ALGEBRAIC[cipaordv1_tau_relp] is tau_relp in component ryr (millisecond).
 * STATES[cipaordv1_Jrelnp] is Jrelnp in component ryr (dimensionless).
 * STATES[cipaordv1_Jrelp] is Jrelp in component ryr (dimensionless).
 * CONSTANTS[cipaordv1_btp] is btp in component ryr (millisecond).
 * CONSTANTS[cipaordv1_a_relp] is a_relp in component ryr (millisecond).
 * ALGEBRAIC[cipaordv1_Jrel_inf_temp] is Jrel_inf_temp in component ryr (dimensionless).
 * ALGEBRAIC[cipaordv1_fJrelp] is fJrelp in component ryr (dimensionless).
 * CONSTANTS[cipaordv1_Jrel_scaling_factor] is Jrel_scaling_factor in component ryr (dimensionless).
 * ALGEBRAIC[cipaordv1_tau_rel_temp] is tau_rel_temp in component ryr (millisecond).
 * ALGEBRAIC[cipaordv1_tau_relp_temp] is tau_relp_temp in component ryr (millisecond).
 * CONSTANTS[cipaordv1_upScale] is upScale in component SERCA (dimensionless).
 * ALGEBRAIC[cipaordv1_Jupnp] is Jupnp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_Jupp] is Jupp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[cipaordv1_fJupp] is fJupp in component SERCA (dimensionless).
 * ALGEBRAIC[cipaordv1_Jleak] is Jleak in component SERCA (millimolar_per_millisecond).
 * CONSTANTS[cipaordv1_Jup_b] is Jup_b in component SERCA (dimensionless).
 * RATES[cipaordv1_V] is d/dt v in component membrane (millivolt).
 * RATES[cipaordv1_CaMKt] is d/dt CaMKt in component CaMK (millimolar).
 * RATES[cipaordv1_nai] is d/dt nai in component intracellular_ions (millimolar).
 * RATES[cipaordv1_nass] is d/dt nass in component intracellular_ions (millimolar).
 * RATES[cipaordv1_ki] is d/dt ki in component intracellular_ions (millimolar).
 * RATES[cipaordv1_kss] is d/dt kss in component intracellular_ions (millimolar).
 * RATES[cipaordv1_cai] is d/dt cai in component intracellular_ions (millimolar).
 * RATES[cipaordv1_cass] is d/dt cass in component intracellular_ions (millimolar).
 * RATES[cipaordv1_cansr] is d/dt cansr in component intracellular_ions (millimolar).
 * RATES[cipaordv1_cajsr] is d/dt cajsr in component intracellular_ions (millimolar).
 * RATES[cipaordv1_m] is d/dt m in component INa (dimensionless).
 * RATES[cipaordv1_hf] is d/dt hf in component INa (dimensionless).
 * RATES[cipaordv1_hs] is d/dt hs in component INa (dimensionless).
 * RATES[cipaordv1_j] is d/dt j in component INa (dimensionless).
 * RATES[cipaordv1_hsp] is d/dt hsp in component INa (dimensionless).
 * RATES[cipaordv1_jp] is d/dt jp in component INa (dimensionless).
 * RATES[cipaordv1_mL] is d/dt mL in component INaL (dimensionless).
 * RATES[cipaordv1_hL] is d/dt hL in component INaL (dimensionless).
 * RATES[cipaordv1_hLp] is d/dt hLp in component INaL (dimensionless).
 * RATES[cipaordv1_a] is d/dt a in component Ito (dimensionless).
 * RATES[cipaordv1_iF] is d/dt iF in component Ito (dimensionless).
 * RATES[cipaordv1_iS] is d/dt iS in component Ito (dimensionless).
 * RATES[cipaordv1_ap] is d/dt ap in component Ito (dimensionless).
 * RATES[cipaordv1_iFp] is d/dt iFp in component Ito (dimensionless).
 * RATES[cipaordv1_iSp] is d/dt iSp in component Ito (dimensionless).
 * RATES[cipaordv1_d] is d/dt d in component ICaL (dimensionless).
 * RATES[cipaordv1_ff] is d/dt ff in component ICaL (dimensionless).
 * RATES[cipaordv1_fs] is d/dt fs in component ICaL (dimensionless).
 * RATES[cipaordv1_fcaf] is d/dt fcaf in component ICaL (dimensionless).
 * RATES[cipaordv1_fcas] is d/dt fcas in component ICaL (dimensionless).
 * RATES[cipaordv1_jca] is d/dt jca in component ICaL (dimensionless).
 * RATES[cipaordv1_ffp] is d/dt ffp in component ICaL (dimensionless).
 * RATES[cipaordv1_fcafp] is d/dt fcafp in component ICaL (dimensionless).
 * RATES[cipaordv1_nca] is d/dt nca in component ICaL (dimensionless).
 * RATES[cipaordv1_IC1] is d/dt IC1 in component IKr (dimensionless).
 * RATES[cipaordv1_IC2] is d/dt IC2 in component IKr (dimensionless).
 * RATES[cipaordv1_C1] is d/dt C1 in component IKr (dimensionless).
 * RATES[cipaordv1_C2] is d/dt C2 in component IKr (dimensionless).
 * RATES[cipaordv1_O] is d/dt O in component IKr (dimensionless).
 * RATES[cipaordv1_IO] is d/dt IO in component IKr (dimensionless).
 * RATES[cipaordv1_IObound] is d/dt IObound in component IKr (dimensionless).
 * RATES[cipaordv1_Obound] is d/dt Obound in component IKr (dimensionless).
 * RATES[cipaordv1_Cbound] is d/dt Cbound in component IKr (dimensionless).
 * RATES[cipaordv1_D] is d/dt D in component IKr (dimensionless).
 * RATES[cipaordv1_xs1] is d/dt xs1 in component IKs (dimensionless).
 * RATES[cipaordv1_xs2] is d/dt xs2 in component IKs (dimensionless).
 * RATES[cipaordv1_xk1] is d/dt xk1 in component IK1 (dimensionless).
 * RATES[cipaordv1_Jrelnp] is d/dt Jrelnp in component ryr (dimensionless).
 * RATES[cipaordv1_Jrelp] is d/dt Jrelp in component ryr (dimensionless).
 */


cipaordv1_land2016::cipaordv1_land2016()
{
algebraic_size = 200 + 24;
constants_size = 206 + 30 + 7 + 1;
states_size = 50 + 7;
ALGEBRAIC = new double[algebraic_size];
CONSTANTS = new double[constants_size];
RATES = new double[states_size];
STATES = new double[states_size];
}

cipaordv1_land2016::~cipaordv1_land2016()
{
delete []ALGEBRAIC;
delete []CONSTANTS;
delete []RATES;
delete []STATES;
}

void cipaordv1_land2016::initConsts()
{
  ___initConsts(0.);
}

void cipaordv1_land2016::initConsts(double type)
{
  ___initConsts(type);
}

void cipaordv1_land2016::initConsts(double type, double conc, const double *hill, const double *herg, bool is_skinned, bool BETA)
{
  std::cout << "initConsts for cipaordv1_land2016\n";
  ___initConsts(type);
  std::cout << "Celltype: " << CONSTANTS[cipaordv1_celltype] << std::endl;
  std::cout << "Control Hill inhibition " << CONSTANTS[cipaordv1_PCa] << " " << CONSTANTS[cipaordv1_GK1] << " " << CONSTANTS[cipaordv1_GKs] << " " << CONSTANTS[cipaordv1_GNaL] << " " << CONSTANTS[cipaordv1_GNa] << " " << CONSTANTS[cipaordv1_Gto] <<std::endl;
  ___applyDrugEffect(conc, hill);
  std::cout << "After drug " << CONSTANTS[cipaordv1_PCa] << " " << CONSTANTS[cipaordv1_GK1] << " " << CONSTANTS[cipaordv1_GKs] << " " << CONSTANTS[cipaordv1_GNaL] << " " << CONSTANTS[cipaordv1_GNa] << " " << CONSTANTS[cipaordv1_Gto] << std::endl;
  std::cout << "Control hERG binding " << CONSTANTS[cipaordv1_Kmax] << " " << CONSTANTS[cipaordv1_Ku] << " " << CONSTANTS[cipaordv1_n] << " " << CONSTANTS[cipaordv1_halfmax] << " " << CONSTANTS[cipaordv1_Vhalf] << " " << CONSTANTS[cipaordv1_Kt] << " " << CONSTANTS[cipaordv1_cnc] << " " << std::endl;
  ___applyHERGBinding(conc, herg);
  std::cout << "After drug " << CONSTANTS[cipaordv1_Kmax] << " " << CONSTANTS[cipaordv1_Ku] << " " << CONSTANTS[cipaordv1_n] << " " << CONSTANTS[cipaordv1_halfmax] << " " << CONSTANTS[cipaordv1_Vhalf] << " " << CONSTANTS[cipaordv1_Kt] << " " << CONSTANTS[cipaordv1_cnc] << " " << std::endl;
  ___initConsts(is_skinned,BETA);
}

void cipaordv1_land2016::___initConsts(double type)
{
CONSTANTS[cipaordv1_celltype] = type;
CONSTANTS[cipaordv1_nao] = 140;
CONSTANTS[cipaordv1_cao] = 1.8;
CONSTANTS[cipaordv1_ko] = 5.4;
CONSTANTS[cipaordv1_R] = 8314;
CONSTANTS[cipaordv1_T] = 310;
CONSTANTS[cipaordv1_F] = 96485;
CONSTANTS[cipaordv1_zna] = 1;
CONSTANTS[cipaordv1_zca] = 2;
CONSTANTS[cipaordv1_zk] = 1;
CONSTANTS[cipaordv1_L] = 0.01;
CONSTANTS[cipaordv1_rad] = 0.0011;
STATES[cipaordv1_V] = -88.00190465;
CONSTANTS[cipaordv1_stim_start] = 10;
CONSTANTS[cipaordv1_stim_end] = 100000000000000000;
CONSTANTS[cipaordv1_amp] = -80;
CONSTANTS[cipaordv1_BCL] = 1000;
CONSTANTS[cipaordv1_duration] = 0.5;
CONSTANTS[cipaordv1_KmCaMK] = 0.15;
CONSTANTS[cipaordv1_aCaMK] = 0.05;
CONSTANTS[cipaordv1_bCaMK] = 0.00068;
CONSTANTS[cipaordv1_CaMKo] = 0.05;
CONSTANTS[cipaordv1_KmCaM] = 0.0015;
STATES[cipaordv1_CaMKt] = 0.0125840447;
STATES[cipaordv1_cass] = 8.49e-05;
CONSTANTS[cipaordv1_cmdnmax_b] = 0.05;
CONSTANTS[cipaordv1_kmcmdn] = 0.00238;
CONSTANTS[cipaordv1_trpnmax] = 0.07;
CONSTANTS[cipaordv1_kmtrpn] = 0.0005;
CONSTANTS[cipaordv1_BSRmax] = 0.047;
CONSTANTS[cipaordv1_KmBSR] = 0.00087;
CONSTANTS[cipaordv1_BSLmax] = 1.124;
CONSTANTS[cipaordv1_KmBSL] = 0.0087;
CONSTANTS[cipaordv1_csqnmax] = 10;
CONSTANTS[cipaordv1_kmcsqn] = 0.8;
STATES[cipaordv1_nai] = 7.268004498;
STATES[cipaordv1_nass] = 7.268089977;
STATES[cipaordv1_ki] = 144.6555918;
STATES[cipaordv1_kss] = 144.6555651;
STATES[cipaordv1_cansr] = 1.619574538;
STATES[cipaordv1_cajsr] = 1.571234014;
STATES[cipaordv1_cai] = 8.6e-05;
CONSTANTS[cipaordv1_cm] = 1;
CONSTANTS[cipaordv1_PKNa] = 0.01833;
CONSTANTS[cipaordv1_mssV1] = 39.57;
CONSTANTS[cipaordv1_mssV2] = 9.871;
CONSTANTS[cipaordv1_mtV1] = 11.64;
CONSTANTS[cipaordv1_mtV2] = 34.77;
CONSTANTS[cipaordv1_mtD1] = 6.765;
CONSTANTS[cipaordv1_mtD2] = 8.552;
CONSTANTS[cipaordv1_mtV3] = 77.42;
CONSTANTS[cipaordv1_mtV4] = 5.955;
STATES[cipaordv1_m] = 0.007344121102;
CONSTANTS[cipaordv1_hssV1] = 82.9;
CONSTANTS[cipaordv1_hssV2] = 6.086;
CONSTANTS[cipaordv1_Ahf] = 0.99;
STATES[cipaordv1_hf] = 0.6981071913;
STATES[cipaordv1_hs] = 0.6980895801;
CONSTANTS[cipaordv1_GNa] = 75;
CONSTANTS[cipaordv1_shift_INa_inact] = 0;
STATES[cipaordv1_j] = 0.6979908432;
STATES[cipaordv1_hsp] = 0.4549485525;
STATES[cipaordv1_jp] = 0.6979245865;
STATES[cipaordv1_mL] = 0.0001882617273;
CONSTANTS[cipaordv1_thL] = 200;
STATES[cipaordv1_hL] = 0.5008548855;
STATES[cipaordv1_hLp] = 0.2693065357;
CONSTANTS[cipaordv1_GNaL_b] = 0.019957499999999975;
CONSTANTS[cipaordv1_Gto_b] = 0.02;
STATES[cipaordv1_a] = 0.001001097687;
STATES[cipaordv1_iF] = 0.9995541745;
STATES[cipaordv1_iS] = 0.5865061736;
STATES[cipaordv1_ap] = 0.0005100862934;
STATES[cipaordv1_iFp] = 0.9995541823;
STATES[cipaordv1_iSp] = 0.6393399482;
CONSTANTS[cipaordv1_Kmn] = 0.002;
CONSTANTS[cipaordv1_k2n] = 1000;
CONSTANTS[cipaordv1_PCa_b] = 0.0001007;
STATES[cipaordv1_d] = 2.34e-9;
STATES[cipaordv1_ff] = 0.9999999909;
STATES[cipaordv1_fs] = 0.9102412777;
STATES[cipaordv1_fcaf] = 0.9999999909;
STATES[cipaordv1_fcas] = 0.9998046777;
STATES[cipaordv1_jca] = 0.9999738312;
STATES[cipaordv1_ffp] = 0.9999999909;
STATES[cipaordv1_fcafp] = 0.9999999909;
STATES[cipaordv1_nca] = 0.002749414044;
CONSTANTS[cipaordv1_GKr_b] = 0.04658545454545456;
STATES[cipaordv1_IC1] = 0.999637;
STATES[cipaordv1_IC2] = 6.83208e-05;
STATES[cipaordv1_C1] = 1.80145e-08;
STATES[cipaordv1_C2] = 8.26619e-05;
STATES[cipaordv1_O] = 0.00015551;
STATES[cipaordv1_IO] = 5.67623e-05;
STATES[cipaordv1_IObound] = 0;
STATES[cipaordv1_Obound] = 0;
STATES[cipaordv1_Cbound] = 0;
STATES[cipaordv1_D] = 0;
CONSTANTS[cipaordv1_A1] = 0.0264;
CONSTANTS[cipaordv1_B1] = 4.631E-05;
CONSTANTS[cipaordv1_q1] = 4.843;
CONSTANTS[cipaordv1_A2] = 4.986E-06;
CONSTANTS[cipaordv1_B2] = -0.004226;
CONSTANTS[cipaordv1_q2] = 4.23;
CONSTANTS[cipaordv1_A3] = 0.001214;
CONSTANTS[cipaordv1_B3] = 0.008516;
CONSTANTS[cipaordv1_q3] = 4.962;
CONSTANTS[cipaordv1_A4] = 1.854E-05;
CONSTANTS[cipaordv1_B4] = -0.04641;
CONSTANTS[cipaordv1_q4] = 3.769;
CONSTANTS[cipaordv1_A11] = 0.0007868;
CONSTANTS[cipaordv1_B11] = 1.535E-08;
CONSTANTS[cipaordv1_q11] = 4.942;
CONSTANTS[cipaordv1_A21] = 5.455E-06;
CONSTANTS[cipaordv1_B21] = -0.1688;
CONSTANTS[cipaordv1_q21] = 4.156;
CONSTANTS[cipaordv1_A31] = 0.005509;
CONSTANTS[cipaordv1_B31] = 7.771E-09;
CONSTANTS[cipaordv1_q31] = 4.22;
CONSTANTS[cipaordv1_A41] = 0.001416;
CONSTANTS[cipaordv1_B41] = -0.02877;
CONSTANTS[cipaordv1_q41] = 1.459;
CONSTANTS[cipaordv1_A51] = 0.4492;
CONSTANTS[cipaordv1_B51] = 0.008595;
CONSTANTS[cipaordv1_q51] = 5;
CONSTANTS[cipaordv1_A52] = 0.3181;
CONSTANTS[cipaordv1_B52] = 3.613E-08;
CONSTANTS[cipaordv1_q52] = 4.663;
CONSTANTS[cipaordv1_A53] = 0.149;
CONSTANTS[cipaordv1_B53] = 0.004668;
CONSTANTS[cipaordv1_q53] = 2.412;
CONSTANTS[cipaordv1_A61] = 0.01241;
CONSTANTS[cipaordv1_B61] = 0.1725;
CONSTANTS[cipaordv1_q61] = 5.568;
CONSTANTS[cipaordv1_A62] = 0.3226;
CONSTANTS[cipaordv1_B62] = -0.0006575;
CONSTANTS[cipaordv1_q62] = 5;
CONSTANTS[cipaordv1_A63] = 0.008978;
CONSTANTS[cipaordv1_B63] = -0.02215;
CONSTANTS[cipaordv1_q63] = 5.682;
CONSTANTS[cipaordv1_Kmax] = 0;
CONSTANTS[cipaordv1_Ku] = 0;
CONSTANTS[cipaordv1_n] = 1;
CONSTANTS[cipaordv1_halfmax] = 1;
CONSTANTS[cipaordv1_Kt] = 0;
CONSTANTS[cipaordv1_Vhalf] = 1;
CONSTANTS[cipaordv1_Temp] = 37;
CONSTANTS[cipaordv1_GKs_b] = 0.006358000000000001;
CONSTANTS[cipaordv1_txs1_max] = 817.3;
STATES[cipaordv1_xs1] = 0.2707758025;
STATES[cipaordv1_xs2] = 0.0001928503426;
CONSTANTS[cipaordv1_GK1_b] = 0.3239783999999998;
STATES[cipaordv1_xk1] = 0.9967597594;
CONSTANTS[cipaordv1_kna1] = 15;
CONSTANTS[cipaordv1_kna2] = 5;
CONSTANTS[cipaordv1_kna3] = 88.12;
CONSTANTS[cipaordv1_kasymm] = 12.5;
CONSTANTS[cipaordv1_wna] = 6e4;
CONSTANTS[cipaordv1_wca] = 6e4;
CONSTANTS[cipaordv1_wnaca] = 5e3;
CONSTANTS[cipaordv1_kcaon] = 1.5e6;
CONSTANTS[cipaordv1_kcaoff] = 5e3;
CONSTANTS[cipaordv1_qna] = 0.5224;
CONSTANTS[cipaordv1_qca] = 0.167;
CONSTANTS[cipaordv1_KmCaAct] = 150e-6;
CONSTANTS[cipaordv1_Gncx_b] = 0.0008;
CONSTANTS[cipaordv1_k1p] = 949.5;
CONSTANTS[cipaordv1_k1m] = 182.4;
CONSTANTS[cipaordv1_k2p] = 687.2;
CONSTANTS[cipaordv1_k2m] = 39.4;
CONSTANTS[cipaordv1_k3p] = 1899;
CONSTANTS[cipaordv1_k3m] = 79300;
CONSTANTS[cipaordv1_k4p] = 639;
CONSTANTS[cipaordv1_k4m] = 40;
CONSTANTS[cipaordv1_Knai0] = 9.073;
CONSTANTS[cipaordv1_Knao0] = 27.78;
CONSTANTS[cipaordv1_delta] = -0.155;
CONSTANTS[cipaordv1_Kki] = 0.5;
CONSTANTS[cipaordv1_Kko] = 0.3582;
CONSTANTS[cipaordv1_MgADP] = 0.05;
CONSTANTS[cipaordv1_MgATP] = 9.8;
CONSTANTS[cipaordv1_Kmgatp] = 1.698e-7;
CONSTANTS[cipaordv1_H] = 1e-7;
CONSTANTS[cipaordv1_eP] = 4.2;
CONSTANTS[cipaordv1_Khp] = 1.698e-7;
CONSTANTS[cipaordv1_Knap] = 224;
CONSTANTS[cipaordv1_Kxkur] = 292;
CONSTANTS[cipaordv1_Pnak_b] = 30;
CONSTANTS[cipaordv1_GKb_b] = 0.003;
CONSTANTS[cipaordv1_PNab] = 3.75e-10;
CONSTANTS[cipaordv1_PCab] = 2.5e-8;
CONSTANTS[cipaordv1_GpCa] = 0.0005;
CONSTANTS[cipaordv1_KmCap] = 0.0005;
CONSTANTS[cipaordv1_bt] = 4.75;
STATES[cipaordv1_Jrelnp] = 2.5e-7;
STATES[cipaordv1_Jrelp] = 3.12e-7;
CONSTANTS[cipaordv1_Jrel_scaling_factor] = 1.0;
CONSTANTS[cipaordv1_Jup_b] = 1.0;
CONSTANTS[cipaordv1_frt] = CONSTANTS[cipaordv1_F]/( CONSTANTS[cipaordv1_R]*CONSTANTS[cipaordv1_T]);
CONSTANTS[cipaordv1_cmdnmax] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_cmdnmax_b]*1.30000 : CONSTANTS[cipaordv1_cmdnmax_b]);
CONSTANTS[cipaordv1_Ahs] = 1.00000 - CONSTANTS[cipaordv1_Ahf];
CONSTANTS[cipaordv1_thLp] =  3.00000*CONSTANTS[cipaordv1_thL];
CONSTANTS[cipaordv1_GNaL] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_GNaL_b]*0.600000 : CONSTANTS[cipaordv1_GNaL_b]);
CONSTANTS[cipaordv1_Gto] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_Gto_b]*4.00000 : CONSTANTS[cipaordv1_celltype]==2.00000 ?  CONSTANTS[cipaordv1_Gto_b]*4.00000 : CONSTANTS[cipaordv1_Gto_b]);
CONSTANTS[cipaordv1_Aff] = 0.600000;
CONSTANTS[cipaordv1_PCa] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_PCa_b]*1.20000 : CONSTANTS[cipaordv1_celltype]==2.00000 ?  CONSTANTS[cipaordv1_PCa_b]*2.50000 : CONSTANTS[cipaordv1_PCa_b]);
CONSTANTS[cipaordv1_tjca] = 75.0000;
CONSTANTS[cipaordv1_v0_CaL] = 0.000000;
CONSTANTS[cipaordv1_GKr] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_GKr_b]*1.30000 : CONSTANTS[cipaordv1_celltype]==2.00000 ?  CONSTANTS[cipaordv1_GKr_b]*0.800000 : CONSTANTS[cipaordv1_GKr_b]);
CONSTANTS[cipaordv1_GKs] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_GKs_b]*1.40000 : CONSTANTS[cipaordv1_GKs_b]);
CONSTANTS[cipaordv1_GK1] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_GK1_b]*1.20000 : CONSTANTS[cipaordv1_celltype]==2.00000 ?  CONSTANTS[cipaordv1_GK1_b]*1.30000 : CONSTANTS[cipaordv1_GK1_b]);
CONSTANTS[cipaordv1_vcell] =  1000.00*3.14000*CONSTANTS[cipaordv1_rad]*CONSTANTS[cipaordv1_rad]*CONSTANTS[cipaordv1_L];
CONSTANTS[cipaordv1_GKb] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_GKb_b]*0.600000 : CONSTANTS[cipaordv1_GKb_b]);
CONSTANTS[cipaordv1_v0_Nab] = 0.000000;
CONSTANTS[cipaordv1_v0_Cab] = 0.000000;
CONSTANTS[cipaordv1_a_rel] =  0.500000*CONSTANTS[cipaordv1_bt];
CONSTANTS[cipaordv1_btp] =  1.25000*CONSTANTS[cipaordv1_bt];
CONSTANTS[cipaordv1_upScale] = (CONSTANTS[cipaordv1_celltype]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[cipaordv1_cnc] = 0.000000;
CONSTANTS[cipaordv1_ffrt] =  CONSTANTS[cipaordv1_F]*CONSTANTS[cipaordv1_frt];
CONSTANTS[cipaordv1_Afs] = 1.00000 - CONSTANTS[cipaordv1_Aff];
CONSTANTS[cipaordv1_PCap] =  1.10000*CONSTANTS[cipaordv1_PCa];
CONSTANTS[cipaordv1_PCaNa] =  0.00125000*CONSTANTS[cipaordv1_PCa];
CONSTANTS[cipaordv1_PCaK] =  0.000357400*CONSTANTS[cipaordv1_PCa];
CONSTANTS[cipaordv1_B_1] =  2.00000*CONSTANTS[cipaordv1_frt];
CONSTANTS[cipaordv1_B_2] = CONSTANTS[cipaordv1_frt];
CONSTANTS[cipaordv1_B_3] = CONSTANTS[cipaordv1_frt];
CONSTANTS[cipaordv1_Ageo] =  2.00000*3.14000*CONSTANTS[cipaordv1_rad]*CONSTANTS[cipaordv1_rad]+ 2.00000*3.14000*CONSTANTS[cipaordv1_rad]*CONSTANTS[cipaordv1_L];
CONSTANTS[cipaordv1_B_Nab] = CONSTANTS[cipaordv1_frt];
CONSTANTS[cipaordv1_B_Cab] =  2.00000*CONSTANTS[cipaordv1_frt];
CONSTANTS[cipaordv1_a_relp] =  0.500000*CONSTANTS[cipaordv1_btp];
CONSTANTS[cipaordv1_PCaNap] =  0.00125000*CONSTANTS[cipaordv1_PCap];
CONSTANTS[cipaordv1_PCaKp] =  0.000357400*CONSTANTS[cipaordv1_PCap];
CONSTANTS[cipaordv1_Acap] =  2.00000*CONSTANTS[cipaordv1_Ageo];
CONSTANTS[cipaordv1_vmyo] =  0.680000*CONSTANTS[cipaordv1_vcell];
CONSTANTS[cipaordv1_vnsr] =  0.0552000*CONSTANTS[cipaordv1_vcell];
CONSTANTS[cipaordv1_vjsr] =  0.00480000*CONSTANTS[cipaordv1_vcell];
CONSTANTS[cipaordv1_vss] =  0.0200000*CONSTANTS[cipaordv1_vcell];
CONSTANTS[cipaordv1_h10_i] = CONSTANTS[cipaordv1_kasymm]+1.00000+ (CONSTANTS[cipaordv1_nao]/CONSTANTS[cipaordv1_kna1])*(1.00000+CONSTANTS[cipaordv1_nao]/CONSTANTS[cipaordv1_kna2]);
CONSTANTS[cipaordv1_h11_i] = ( CONSTANTS[cipaordv1_nao]*CONSTANTS[cipaordv1_nao])/( CONSTANTS[cipaordv1_h10_i]*CONSTANTS[cipaordv1_kna1]*CONSTANTS[cipaordv1_kna2]);
CONSTANTS[cipaordv1_h12_i] = 1.00000/CONSTANTS[cipaordv1_h10_i];
CONSTANTS[cipaordv1_k1_i] =  CONSTANTS[cipaordv1_h12_i]*CONSTANTS[cipaordv1_cao]*CONSTANTS[cipaordv1_kcaon];
CONSTANTS[cipaordv1_k2_i] = CONSTANTS[cipaordv1_kcaoff];
CONSTANTS[cipaordv1_k5_i] = CONSTANTS[cipaordv1_kcaoff];
CONSTANTS[cipaordv1_Gncx] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_Gncx_b]*1.10000 : CONSTANTS[cipaordv1_celltype]==2.00000 ?  CONSTANTS[cipaordv1_Gncx_b]*1.40000 : CONSTANTS[cipaordv1_Gncx_b]);
CONSTANTS[cipaordv1_h10_ss] = CONSTANTS[cipaordv1_kasymm]+1.00000+ (CONSTANTS[cipaordv1_nao]/CONSTANTS[cipaordv1_kna1])*(1.00000+CONSTANTS[cipaordv1_nao]/CONSTANTS[cipaordv1_kna2]);
CONSTANTS[cipaordv1_h11_ss] = ( CONSTANTS[cipaordv1_nao]*CONSTANTS[cipaordv1_nao])/( CONSTANTS[cipaordv1_h10_ss]*CONSTANTS[cipaordv1_kna1]*CONSTANTS[cipaordv1_kna2]);
CONSTANTS[cipaordv1_h12_ss] = 1.00000/CONSTANTS[cipaordv1_h10_ss];
CONSTANTS[cipaordv1_k1_ss] =  CONSTANTS[cipaordv1_h12_ss]*CONSTANTS[cipaordv1_cao]*CONSTANTS[cipaordv1_kcaon];
CONSTANTS[cipaordv1_k2_ss] = CONSTANTS[cipaordv1_kcaoff];
CONSTANTS[cipaordv1_k5_ss] = CONSTANTS[cipaordv1_kcaoff];
CONSTANTS[cipaordv1_b1] =  CONSTANTS[cipaordv1_k1m]*CONSTANTS[cipaordv1_MgADP];
CONSTANTS[cipaordv1_a2] = CONSTANTS[cipaordv1_k2p];
CONSTANTS[cipaordv1_a4] = (( CONSTANTS[cipaordv1_k4p]*CONSTANTS[cipaordv1_MgATP])/CONSTANTS[cipaordv1_Kmgatp])/(1.00000+CONSTANTS[cipaordv1_MgATP]/CONSTANTS[cipaordv1_Kmgatp]);
CONSTANTS[cipaordv1_Pnak] = (CONSTANTS[cipaordv1_celltype]==1.00000 ?  CONSTANTS[cipaordv1_Pnak_b]*0.900000 : CONSTANTS[cipaordv1_celltype]==2.00000 ?  CONSTANTS[cipaordv1_Pnak_b]*0.700000 : CONSTANTS[cipaordv1_Pnak_b]);
}

void cipaordv1_land2016::___initConsts(bool is_skinned, bool BETA){
    STATES[land2016_XS] = 0.0;
    STATES[land2016_XW] = 0.0;
    STATES[land2016_TRPN] = 0.0;
    STATES[land2016_TmBlocked] = 0.0;
    STATES[land2016_ZETAS] = 0.0;
    STATES[land2016_ZETAW] = 0.0;
    STATES[land2016_dCd_dt] = 0.0;
	    // user input
    CONSTANTS[land2016_dlambda_dt] = 0; 
    CONSTANTS[land2016_lambda] = 1.0;
    CONSTANTS[land2016_Cai] = 0.0;
    // if (t == -1) {
        RATES[land2016_XS] = 0;
        RATES[land2016_XW] = 0;
        RATES[land2016_TRPN] = 0;
        RATES[land2016_TmBlocked] = 1;
        RATES[land2016_ZETAS] = 0;
        RATES[land2016_ZETAW] = 0;
        RATES[land2016_dCd_dt] = 0;
    //     return;
    // }

	if (CONSTANTS[land2016_lambda] >= 1.2){
		CONSTANTS[land2016_lambda] = 1.2;
	}
    CONSTANTS[land2016_perm50] = 0.35;
    CONSTANTS[land2016_TRPN_n] = 2;
    CONSTANTS[land2016_koff] = 0.1;
    CONSTANTS[land2016_dr] = 0.25;
    CONSTANTS[land2016_wfrac] = 0.5;

    CONSTANTS[land2016_TOT_A] = 25;
    // CONSTANTS[land2016_ktm_unblock] = 1; // before margara
    CONSTANTS[land2016_ktm_unblock] = 0.04; // after margara

    CONSTANTS[land2016_beta_1] = -2.4;
    CONSTANTS[land2016_beta_0] = 2.3;
    CONSTANTS[land2016_gamma_idx] = 0.0085;
    CONSTANTS[land2016_gamma_wu] = 0.615;
    CONSTANTS[land2016_phi] = 2.23;  

    if (is_skinned == true){
      CONSTANTS[land2016_nperm] = 2.2;
      ALGEBRAIC[land2016_ca50] = 2.5;
      CONSTANTS[land2016_Tref] = 40.5;
      CONSTANTS[land2016_nu] = 1;
      CONSTANTS[land2016_mu] = 1;
    } else{
      CONSTANTS[land2016_nperm] = 5.0; // before margara
      CONSTANTS[land2016_nperm] = 2.4; // after margara
      ALGEBRAIC[land2016_ca50] = 0.805;
      CONSTANTS[land2016_Tref] = 120.0;
      CONSTANTS[land2016_nu] = 7;
      CONSTANTS[land2016_mu] = 3; 
    }
    if (BETA == true){
        // // input an array of beta, 
        // CONSTANTS[land2016_beta_1] = beta[1];
        // CONSTANTS[land2016_beta_0] = beta[0];
    }
    
    // k_ws = 0.004 * mu;
    // k_uw = 0.026  * nu;
    CONSTANTS[land2016_k_ws] = 0.004 * CONSTANTS[land2016_mu];
    CONSTANTS[land2016_k_uw] = 0.026 * CONSTANTS[land2016_nu];

    //hyperparameters for land2016
    CONSTANTS[land2016_y1] = 0.0;
    CONSTANTS[land2016_y2] = 0.0;
    CONSTANTS[land2016_y3] = 0.0;
    CONSTANTS[land2016_y4] = 1.0;
    CONSTANTS[land2016_y5] = 0.0;
    CONSTANTS[land2016_y6] = 0.0;
    CONSTANTS[land2016_y7] = 0.0;

    //STATES Variables
    if(CONSTANTS[land2016_y1]>0) STATES[land2016_XS] = fmax(0,CONSTANTS[land2016_y1]);
    if(CONSTANTS[land2016_y2]>0) STATES[land2016_XW] = fmax(0,CONSTANTS[land2016_y2]);
    if(CONSTANTS[land2016_y3]>0) STATES[land2016_TRPN] = fmax(0,CONSTANTS[land2016_y3]);
    STATES[land2016_TmBlocked] = CONSTANTS[land2016_y4];
    STATES[land2016_ZETAS] = CONSTANTS[land2016_y5] ;
    STATES[land2016_ZETAW] = CONSTANTS[land2016_y6] ;
    STATES[land2016_dCd_dt] = CONSTANTS[land2016_y7] ;
    // printf("initialisation\n");


    // seventh state is below, in passive model (dCd_dt)

    CONSTANTS[land2016_par_k] = 7;
    CONSTANTS[land2016_b] = 9.1;
    CONSTANTS[land2016_eta_l] = 200;
    CONSTANTS[land2016_eta_s] = 20;
    CONSTANTS[land2016_land_a] = 2.1;

    // cai and cai scaling
    CONSTANTS[land2016_cai_scaling] = 0;
}

void cipaordv1_land2016::computeRates( double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC ){
// passing the calcium information from cipaordv1 to land2016
if (CONSTANTS[land2016_cai_scaling] == 1){
  CONSTANTS[land2016_Cai] = STATES[cipaordv1_cai]*1000.;
} else{
  CONSTANTS[land2016_Cai] = STATES[cipaordv1_cai];
}
//=========
// land2016
//=========

// XB model

// lambda = min(1.2,lambda);
// Lfac =  max(0, 1 + beta_0 * (lambda + min(0.87,lambda) - 1.87) );

CONSTANTS[land2016_lambda] = fmin(1.2, CONSTANTS[land2016_lambda]);
ALGEBRAIC[land2016_Lfac]  = fmax(0, 1 + CONSTANTS[land2016_beta_0] * (CONSTANTS[land2016_lambda] + fmin(0.87, CONSTANTS[land2016_lambda]) - 1.87));

if (ALGEBRAIC[land2016_Lfac] < 0){
    ALGEBRAIC[land2016_Lfac] = 0 ;
}
// cdw = phi * k_uw * (1-dr)*(1-wfrac) /  ((1-dr)*wfrac);
// cds = phi * k_ws * (1-dr)*wfrac / dr;
ALGEBRAIC[land2016_cdw] = CONSTANTS[land2016_phi] * CONSTANTS[land2016_k_uw] * (1-CONSTANTS[land2016_dr])*(1-CONSTANTS[land2016_wfrac]) /  ((1-CONSTANTS[land2016_dr])*CONSTANTS[land2016_wfrac]);
ALGEBRAIC[land2016_cds] = CONSTANTS[land2016_phi] * CONSTANTS[land2016_k_ws] * (1-CONSTANTS[land2016_dr])*CONSTANTS[land2016_wfrac] / CONSTANTS[land2016_dr];

// k_wu = k_uw * (1/wfrac - 1) - k_ws;
// k_su = k_ws * (1/dr - 1) * wfrac; 
// A = (0.25 * TOT_A) / ((1-dr)*wfrac + dr) * (dr/0.25);
ALGEBRAIC[land2016_k_wu] = CONSTANTS[land2016_k_uw] * (1/CONSTANTS[land2016_wfrac] - 1) - CONSTANTS[land2016_k_ws];
ALGEBRAIC[land2016_k_su] = CONSTANTS[land2016_k_ws] * (1/CONSTANTS[land2016_dr] - 1) * CONSTANTS[land2016_wfrac]; 
ALGEBRAIC[land2016_A] = (0.25 * CONSTANTS[land2016_TOT_A]) / ((1-CONSTANTS[land2016_dr])*CONSTANTS[land2016_wfrac] + CONSTANTS[land2016_dr]) * (CONSTANTS[land2016_dr]/0.25);

// XU   = (1 - TmBlocked) - XW - XS;
ALGEBRAIC[land2016_XU] = (1 - STATES[land2016_TmBlocked]) - STATES[land2016_XW] - STATES[land2016_XS]; 
// unattached available xb = all - tm blocked - already prepowerstroke - already post-poststroke - no overlap

// xb_ws = k_ws * XW;
// xb_uw = k_uw * XU ;
// xb_wu = k_wu * XW;
// xb_su = k_su * XS;
ALGEBRAIC[land2016_xb_ws] = CONSTANTS[land2016_k_ws] * STATES[land2016_XW];
ALGEBRAIC[land2016_xb_uw] = CONSTANTS[land2016_k_uw] * ALGEBRAIC[land2016_XU];
ALGEBRAIC[land2016_xb_wu] = ALGEBRAIC[land2016_k_wu] * STATES[land2016_XW];
ALGEBRAIC[land2016_xb_su] = ALGEBRAIC[land2016_k_su] * STATES[land2016_XS];

//gamma_rate  = gamma * max( (ZETAS > 0) .* ZETAS , (ZETAS < -1) .* (-ZETAS-1) );
// element wise multiplication:
double temp_zetas1, temp_zetas2;
temp_zetas1=temp_zetas2=0;
if(STATES[land2016_ZETAS] > 0){
    temp_zetas1 = STATES[land2016_ZETAS];
}
if (STATES[land2016_ZETAS]<-1){
    temp_zetas2 = -1 * STATES[land2016_ZETAS] -1; 
}
ALGEBRAIC[land2016_gamma_rate] = CONSTANTS[land2016_gamma_idx] * fmax( temp_zetas1, temp_zetas2);

// xb_su_gamma = gamma_rate * XS;
ALGEBRAIC[land2016_xb_su_gamma] = ALGEBRAIC[land2016_gamma_rate] * STATES[land2016_XS];

// gamma_rate_w  = gamma_wu * abs(ZETAW);
ALGEBRAIC[land2016_gamma_rate_w] = CONSTANTS[land2016_gamma_wu] * fabs(STATES[land2016_ZETAW]); //weak xbs don't like being strained
// xb_wu_gamma = gamma_rate_w * XW; 
ALGEBRAIC[land2016_xb_wu_gamma] = ALGEBRAIC[land2016_gamma_rate_w] * STATES[land2016_XW];

// dydt(1)  = xb_ws - xb_su - xb_su_gamma;
// dydt(2)  = xb_uw - xb_wu - xb_ws - xb_wu_gamma;
RATES[land2016_XS] = ALGEBRAIC[land2016_xb_ws] - ALGEBRAIC[land2016_xb_su] - ALGEBRAIC[land2016_xb_wu_gamma];
RATES[land2016_XW] = ALGEBRAIC[land2016_xb_uw] - ALGEBRAIC[land2016_xb_wu] - ALGEBRAIC[land2016_xb_ws] - ALGEBRAIC[land2016_xb_wu_gamma];

// ca50 = ca50 + beta_1*min(0.2,lambda - 1);
// dydt(3)  = koff * ( (Cai/ca50)^TRPN_n * (1-TRPN) - TRPN);
ALGEBRAIC[land2016_ca50] = ALGEBRAIC[land2016_ca50] + CONSTANTS[land2016_beta_1] * fmin(0.2, CONSTANTS[land2016_lambda] - 1);

RATES[land2016_TRPN] = CONSTANTS[land2016_koff] * (pow((CONSTANTS[land2016_Cai] / ALGEBRAIC[land2016_ca50]), CONSTANTS[land2016_TRPN_n]) * (1 - STATES[land2016_TRPN]) - STATES[land2016_TRPN]);
////// Cai constants will trigger NaN, followed by states TRPN

// XSSS = dr * 0.5;
// XWSS = (1-dr) * wfrac * 0.5;
// ktm_block = ktm_unblock * (perm50^nperm) *  0.5 / (0.5 - XSSS - XWSS);
ALGEBRAIC[land2016_XSSS] = CONSTANTS[land2016_dr] * 0.5;
ALGEBRAIC[land2016_XWSS] = (1 - CONSTANTS[land2016_dr]) * CONSTANTS[land2016_wfrac] * 0.5;
ALGEBRAIC[land2016_ktm_block]= CONSTANTS[land2016_ktm_unblock] * (pow(CONSTANTS[land2016_perm50], CONSTANTS[land2016_nperm]) * 0.5) / (0.5 - ALGEBRAIC[land2016_XSSS] - ALGEBRAIC[land2016_XWSS]);

// dydt(4)  = ktm_block * min(100, (TRPN^-(nperm/2))) * XU  - ktm_unblock * (TRPN^(nperm/2)) *  TmBlocked;
RATES[land2016_TmBlocked] = CONSTANTS[land2016_ktm_block] * fmin(100, pow(STATES[land2016_TRPN], -(CONSTANTS[land2016_nperm] / 2))) * ALGEBRAIC[land2016_XU] - CONSTANTS[land2016_ktm_unblock]  * pow(STATES[land2016_TRPN], (CONSTANTS[land2016_nperm] / 2)) * STATES[land2016_TmBlocked];
//-------------------------------------------------------------------------------
// Velocity dependence -- assumes distortion resets on W->S
// dydt(5) = A * dlambda_dt - cds * ZETAS;% - gamma_rate * ZETAS;
// dydt(6) = A * dlambda_dt - cdw * ZETAW;% - gamma_rate_w * ZETAW;
RATES[land2016_ZETAS] = CONSTANTS[land2016_A] * CONSTANTS[land2016_dlambda_dt] - ALGEBRAIC[land2016_cds] * STATES[land2016_ZETAS];
RATES[land2016_ZETAW] = CONSTANTS[land2016_A] * CONSTANTS[land2016_dlambda_dt] - ALGEBRAIC[land2016_cdw] * STATES[land2016_ZETAW];

//-------------------------------------------------------------------------------
// Passive model
// this is quite scary, did not fix the problem but Cd actually y[6] aka y(7)
//-------------------------------------------------------------------------------

// Cd = y(7);
// C = lambda - 1;

// CONSTANTS[land2016_Cd] = STATES[land2016_dCd_dt]; 
CONSTANTS[land2016_Cd] = CONSTANTS[land2016_y7];
CONSTANTS[land2016_C] = CONSTANTS[land2016_lambda] - 1;

// if (C - Cd) < 0
//  eta = eta_s;
// else
//  eta = eta_l;
// end
if ((CONSTANTS[land2016_C] - CONSTANTS[land2016_Cd]) < 0) {
	CONSTANTS[land2016_eta] = CONSTANTS[land2016_eta_s];
} else {
	CONSTANTS[land2016_eta] = CONSTANTS[land2016_eta_l];
}

// dCd_dt = par_k * (C - Cd) / eta; % F2=Fd
STATES[land2016_dCd_dt] = CONSTANTS[land2016_par_k] * (CONSTANTS[land2016_C] - CONSTANTS[land2016_Cd]) / CONSTANTS[land2016_eta];
// dydt(7) = dCd_dt;
RATES[land2016_dCd_dt] = STATES[land2016_dCd_dt];

// Fd = eta * dCd_dt;
ALGEBRAIC[land2016_Fd] = CONSTANTS[land2016_eta] * STATES[land2016_dCd_dt];
// F1 = (exp( b * C) - 1);
ALGEBRAIC[land2016_F1] = (exp(CONSTANTS[land2016_b] * CONSTANTS[land2016_C]) - 1);
// Tp = a * (F1 + Fd);
ALGEBRAIC[land2016_Tp] = CONSTANTS[land2016_land_a] * (ALGEBRAIC[land2016_F1] + ALGEBRAIC[land2016_Fd]);

//-------------------------------------------------------------------------------
// Active and Total Force
//-------------------------------------------------------------------------------

// Ta = Lfac * (Tref/dr) * ( (ZETAS+1) * XS + (ZETAW) * XW );
// T = Ta + Tp;
ALGEBRAIC[land2016_Ta] = ALGEBRAIC[land2016_Lfac] * (CONSTANTS[land2016_Tref] / CONSTANTS[land2016_dr]) * ( (STATES[land2016_ZETAS] + 1) * STATES[land2016_XS] + (STATES[land2016_ZETAW]) * STATES[land2016_XW]);
ALGEBRAIC[land2016_land_T] = ALGEBRAIC[land2016_Ta] + ALGEBRAIC[land2016_Tp];

//=========
// cipaordv1
//=========
// compute all rates
ALGEBRAIC[cipaordv1_hLss] = 1.00000/(1.00000+exp((STATES[cipaordv1_V]+87.6100)/7.48800));
ALGEBRAIC[cipaordv1_hLssp] = 1.00000/(1.00000+exp((STATES[cipaordv1_V]+93.8100)/7.48800));
ALGEBRAIC[cipaordv1_mss] = 1.00000/(1.00000+exp(- (STATES[cipaordv1_V]+CONSTANTS[cipaordv1_mssV1])/CONSTANTS[cipaordv1_mssV2]));
ALGEBRAIC[cipaordv1_tm] = 1.00000/( CONSTANTS[cipaordv1_mtD1]*exp((STATES[cipaordv1_V]+CONSTANTS[cipaordv1_mtV1])/CONSTANTS[cipaordv1_mtV2])+ CONSTANTS[cipaordv1_mtD2]*exp(- (STATES[cipaordv1_V]+CONSTANTS[cipaordv1_mtV3])/CONSTANTS[cipaordv1_mtV4]));
ALGEBRAIC[cipaordv1_hss] = 1.00000/(1.00000+exp(((STATES[cipaordv1_V]+CONSTANTS[cipaordv1_hssV1]) - CONSTANTS[cipaordv1_shift_INa_inact])/CONSTANTS[cipaordv1_hssV2]));
ALGEBRAIC[cipaordv1_thf] = 1.00000/( 1.43200e-05*exp(- ((STATES[cipaordv1_V]+1.19600) - CONSTANTS[cipaordv1_shift_INa_inact])/6.28500)+ 6.14900*exp(((STATES[cipaordv1_V]+0.509600) - CONSTANTS[cipaordv1_shift_INa_inact])/20.2700));
ALGEBRAIC[cipaordv1_ths] = 1.00000/( 0.00979400*exp(- ((STATES[cipaordv1_V]+17.9500) - CONSTANTS[cipaordv1_shift_INa_inact])/28.0500)+ 0.334300*exp(((STATES[cipaordv1_V]+5.73000) - CONSTANTS[cipaordv1_shift_INa_inact])/56.6600));
ALGEBRAIC[cipaordv1_ass] = 1.00000/(1.00000+exp(- (STATES[cipaordv1_V] - 14.3400)/14.8200));
ALGEBRAIC[cipaordv1_ta] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- (STATES[cipaordv1_V] - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[cipaordv1_V]+100.000)/29.3814)));
ALGEBRAIC[cipaordv1_dss] = 1.00000/(1.00000+exp(- (STATES[cipaordv1_V]+3.94000)/4.23000));
ALGEBRAIC[cipaordv1_td] = 0.600000+1.00000/(exp( - 0.0500000*(STATES[cipaordv1_V]+6.00000))+exp( 0.0900000*(STATES[cipaordv1_V]+14.0000)));
ALGEBRAIC[cipaordv1_fss] = 1.00000/(1.00000+exp((STATES[cipaordv1_V]+19.5800)/3.69600));
ALGEBRAIC[cipaordv1_tff] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[cipaordv1_V]+20.0000)/10.0000)+ 0.00450000*exp((STATES[cipaordv1_V]+20.0000)/10.0000));
ALGEBRAIC[cipaordv1_tfs] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[cipaordv1_V]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[cipaordv1_V]+5.00000)/6.00000));
ALGEBRAIC[cipaordv1_fcass] = ALGEBRAIC[cipaordv1_fss];
ALGEBRAIC[cipaordv1_km2n] =  STATES[cipaordv1_jca]*1.00000;
ALGEBRAIC[cipaordv1_anca] = 1.00000/(CONSTANTS[cipaordv1_k2n]/ALGEBRAIC[cipaordv1_km2n]+pow(1.00000+CONSTANTS[cipaordv1_Kmn]/STATES[cipaordv1_cass], 4.00000));
ALGEBRAIC[cipaordv1_xs1ss] = 1.00000/(1.00000+exp(- (STATES[cipaordv1_V]+11.6000)/8.93200));
ALGEBRAIC[cipaordv1_txs1] = CONSTANTS[cipaordv1_txs1_max]+1.00000/( 0.000232600*exp((STATES[cipaordv1_V]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[cipaordv1_V]+210.000)/230.000));
ALGEBRAIC[cipaordv1_xk1ss] = 1.00000/(1.00000+exp(- (STATES[cipaordv1_V]+ 2.55380*CONSTANTS[cipaordv1_ko]+144.590)/( 1.56920*CONSTANTS[cipaordv1_ko]+3.81150)));
ALGEBRAIC[cipaordv1_txk1] = 122.200/(exp(- (STATES[cipaordv1_V]+127.200)/20.3600)+exp((STATES[cipaordv1_V]+236.800)/69.3300));
ALGEBRAIC[cipaordv1_CaMKb] = ( CONSTANTS[cipaordv1_CaMKo]*(1.00000 - STATES[cipaordv1_CaMKt]))/(1.00000+CONSTANTS[cipaordv1_KmCaM]/STATES[cipaordv1_cass]);
ALGEBRAIC[cipaordv1_jss] = ALGEBRAIC[cipaordv1_hss];
ALGEBRAIC[cipaordv1_tj] = 2.03800+1.00000/( 0.0213600*exp(- ((STATES[cipaordv1_V]+100.600) - CONSTANTS[cipaordv1_shift_INa_inact])/8.28100)+ 0.305200*exp(((STATES[cipaordv1_V]+0.994100) - CONSTANTS[cipaordv1_shift_INa_inact])/38.4500));
ALGEBRAIC[cipaordv1_assp] = 1.00000/(1.00000+exp(- (STATES[cipaordv1_V] - 24.3400)/14.8200));
ALGEBRAIC[cipaordv1_tfcaf] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[cipaordv1_V] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[cipaordv1_V] - 4.00000)/7.00000));
ALGEBRAIC[cipaordv1_tfcas] = 100.000+1.00000/( 0.000120000*exp(- STATES[cipaordv1_V]/3.00000)+ 0.000120000*exp(STATES[cipaordv1_V]/7.00000));
ALGEBRAIC[cipaordv1_tffp] =  2.50000*ALGEBRAIC[cipaordv1_tff];
ALGEBRAIC[cipaordv1_xs2ss] = ALGEBRAIC[cipaordv1_xs1ss];
ALGEBRAIC[cipaordv1_txs2] = 1.00000/( 0.0100000*exp((STATES[cipaordv1_V] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[cipaordv1_V]+66.5400)/31.0000));
ALGEBRAIC[cipaordv1_hssp] = 1.00000/(1.00000+exp(((STATES[cipaordv1_V]+89.1000) - CONSTANTS[cipaordv1_shift_INa_inact])/6.08600));
ALGEBRAIC[cipaordv1_thsp] =  3.00000*ALGEBRAIC[cipaordv1_ths];
ALGEBRAIC[cipaordv1_tjp] =  1.46000*ALGEBRAIC[cipaordv1_tj];
ALGEBRAIC[cipaordv1_mLss] = 1.00000/(1.00000+exp(- (STATES[cipaordv1_V]+42.8500)/5.26400));
ALGEBRAIC[cipaordv1_tmL] = ALGEBRAIC[cipaordv1_tm];
ALGEBRAIC[cipaordv1_tfcafp] =  2.50000*ALGEBRAIC[cipaordv1_tfcaf];
ALGEBRAIC[cipaordv1_iss] = 1.00000/(1.00000+exp((STATES[cipaordv1_V]+43.9400)/5.71100));
ALGEBRAIC[cipaordv1_delta_epi] = (CONSTANTS[cipaordv1_celltype]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[cipaordv1_V]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[cipaordv1_tiF_b] = 4.56200+1.00000/( 0.393300*exp(- (STATES[cipaordv1_V]+100.000)/100.000)+ 0.0800400*exp((STATES[cipaordv1_V]+50.0000)/16.5900));
ALGEBRAIC[cipaordv1_tiF] =  ALGEBRAIC[cipaordv1_tiF_b]*ALGEBRAIC[cipaordv1_delta_epi];
ALGEBRAIC[cipaordv1_tiS_b] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[cipaordv1_V]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[cipaordv1_V]+114.100)/8.07900));
ALGEBRAIC[cipaordv1_tiS] =  ALGEBRAIC[cipaordv1_tiS_b]*ALGEBRAIC[cipaordv1_delta_epi];
ALGEBRAIC[cipaordv1_dti_develop] = 1.35400+0.000100000/(exp((STATES[cipaordv1_V] - 167.400)/15.8900)+exp(- (STATES[cipaordv1_V] - 12.2300)/0.215400));
ALGEBRAIC[cipaordv1_dti_recover] = 1.00000 - 0.500000/(1.00000+exp((STATES[cipaordv1_V]+70.0000)/20.0000));
ALGEBRAIC[cipaordv1_tiFp] =  ALGEBRAIC[cipaordv1_dti_develop]*ALGEBRAIC[cipaordv1_dti_recover]*ALGEBRAIC[cipaordv1_tiF];
ALGEBRAIC[cipaordv1_tiSp] =  ALGEBRAIC[cipaordv1_dti_develop]*ALGEBRAIC[cipaordv1_dti_recover]*ALGEBRAIC[cipaordv1_tiS];
ALGEBRAIC[cipaordv1_f] =  CONSTANTS[cipaordv1_Aff]*STATES[cipaordv1_ff]+ CONSTANTS[cipaordv1_Afs]*STATES[cipaordv1_fs];
ALGEBRAIC[cipaordv1_Afcaf] = 0.300000+0.600000/(1.00000+exp((STATES[cipaordv1_V] - 10.0000)/10.0000));
ALGEBRAIC[cipaordv1_Afcas] = 1.00000 - ALGEBRAIC[cipaordv1_Afcaf];
ALGEBRAIC[cipaordv1_fca] =  ALGEBRAIC[cipaordv1_Afcaf]*STATES[cipaordv1_fcaf]+ ALGEBRAIC[cipaordv1_Afcas]*STATES[cipaordv1_fcas];
ALGEBRAIC[cipaordv1_fp] =  CONSTANTS[cipaordv1_Aff]*STATES[cipaordv1_ffp]+ CONSTANTS[cipaordv1_Afs]*STATES[cipaordv1_fs];
ALGEBRAIC[cipaordv1_fcap] =  ALGEBRAIC[cipaordv1_Afcaf]*STATES[cipaordv1_fcafp]+ ALGEBRAIC[cipaordv1_Afcas]*STATES[cipaordv1_fcas];
ALGEBRAIC[cipaordv1_vfrt] =  STATES[cipaordv1_V]*CONSTANTS[cipaordv1_frt];
ALGEBRAIC[cipaordv1_A_1] = ( 4.00000*CONSTANTS[cipaordv1_ffrt]*( STATES[cipaordv1_cass]*exp( 2.00000*ALGEBRAIC[cipaordv1_vfrt]) -  0.341000*CONSTANTS[cipaordv1_cao]))/CONSTANTS[cipaordv1_B_1];
ALGEBRAIC[cipaordv1_U_1] =  CONSTANTS[cipaordv1_B_1]*(STATES[cipaordv1_V] - CONSTANTS[cipaordv1_v0_CaL]);
ALGEBRAIC[cipaordv1_PhiCaL] = (- 1.00000e-07<=ALGEBRAIC[cipaordv1_U_1]&&ALGEBRAIC[cipaordv1_U_1]<=1.00000e-07 ?  ALGEBRAIC[cipaordv1_A_1]*(1.00000 -  0.500000*ALGEBRAIC[cipaordv1_U_1]) : ( ALGEBRAIC[cipaordv1_A_1]*ALGEBRAIC[cipaordv1_U_1])/(exp(ALGEBRAIC[cipaordv1_U_1]) - 1.00000));
ALGEBRAIC[cipaordv1_CaMKa] = ALGEBRAIC[cipaordv1_CaMKb]+STATES[cipaordv1_CaMKt];
ALGEBRAIC[cipaordv1_fICaLp] = 1.00000/(1.00000+CONSTANTS[cipaordv1_KmCaMK]/ALGEBRAIC[cipaordv1_CaMKa]);
ALGEBRAIC[cipaordv1_ICaL] =  (1.00000 - ALGEBRAIC[cipaordv1_fICaLp])*CONSTANTS[cipaordv1_PCa]*ALGEBRAIC[cipaordv1_PhiCaL]*STATES[cipaordv1_d]*( ALGEBRAIC[cipaordv1_f]*(1.00000 - STATES[cipaordv1_nca])+ STATES[cipaordv1_jca]*ALGEBRAIC[cipaordv1_fca]*STATES[cipaordv1_nca])+ ALGEBRAIC[cipaordv1_fICaLp]*CONSTANTS[cipaordv1_PCap]*ALGEBRAIC[cipaordv1_PhiCaL]*STATES[cipaordv1_d]*( ALGEBRAIC[cipaordv1_fp]*(1.00000 - STATES[cipaordv1_nca])+ STATES[cipaordv1_jca]*ALGEBRAIC[cipaordv1_fcap]*STATES[cipaordv1_nca]);
ALGEBRAIC[cipaordv1_Jrel_inf_temp] = ( CONSTANTS[cipaordv1_a_rel]*- ALGEBRAIC[cipaordv1_ICaL])/(1.00000+ 1.00000*pow(1.50000/STATES[cipaordv1_cajsr], 8.00000));
ALGEBRAIC[cipaordv1_Jrel_inf] = (CONSTANTS[cipaordv1_celltype]==2.00000 ?  ALGEBRAIC[cipaordv1_Jrel_inf_temp]*1.70000 : ALGEBRAIC[cipaordv1_Jrel_inf_temp]);
ALGEBRAIC[cipaordv1_tau_rel_temp] = CONSTANTS[cipaordv1_bt]/(1.00000+0.0123000/STATES[cipaordv1_cajsr]);
ALGEBRAIC[cipaordv1_tau_rel] = (ALGEBRAIC[cipaordv1_tau_rel_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[cipaordv1_tau_rel_temp]);
ALGEBRAIC[cipaordv1_Jrel_temp] = ( CONSTANTS[cipaordv1_a_relp]*- ALGEBRAIC[cipaordv1_ICaL])/(1.00000+pow(1.50000/STATES[cipaordv1_cajsr], 8.00000));
ALGEBRAIC[cipaordv1_Jrel_infp] = (CONSTANTS[cipaordv1_celltype]==2.00000 ?  ALGEBRAIC[cipaordv1_Jrel_temp]*1.70000 : ALGEBRAIC[cipaordv1_Jrel_temp]);
ALGEBRAIC[cipaordv1_tau_relp_temp] = CONSTANTS[cipaordv1_btp]/(1.00000+0.0123000/STATES[cipaordv1_cajsr]);
ALGEBRAIC[cipaordv1_tau_relp] = (ALGEBRAIC[cipaordv1_tau_relp_temp]<0.00100000 ? 0.00100000 : ALGEBRAIC[cipaordv1_tau_relp_temp]);
ALGEBRAIC[cipaordv1_EK] =  (( CONSTANTS[cipaordv1_R]*CONSTANTS[cipaordv1_T])/CONSTANTS[cipaordv1_F])*log(CONSTANTS[cipaordv1_ko]/STATES[cipaordv1_ki]);
ALGEBRAIC[cipaordv1_AiF] = 1.00000/(1.00000+exp((STATES[cipaordv1_V] - 213.600)/151.200));
ALGEBRAIC[cipaordv1_AiS] = 1.00000 - ALGEBRAIC[cipaordv1_AiF];
ALGEBRAIC[cipaordv1_i] =  ALGEBRAIC[cipaordv1_AiF]*STATES[cipaordv1_iF]+ ALGEBRAIC[cipaordv1_AiS]*STATES[cipaordv1_iS];
ALGEBRAIC[cipaordv1_ip] =  ALGEBRAIC[cipaordv1_AiF]*STATES[cipaordv1_iFp]+ ALGEBRAIC[cipaordv1_AiS]*STATES[cipaordv1_iSp];
ALGEBRAIC[cipaordv1_fItop] = 1.00000/(1.00000+CONSTANTS[cipaordv1_KmCaMK]/ALGEBRAIC[cipaordv1_CaMKa]);
ALGEBRAIC[cipaordv1_Ito] =  CONSTANTS[cipaordv1_Gto]*(STATES[cipaordv1_V] - ALGEBRAIC[cipaordv1_EK])*( (1.00000 - ALGEBRAIC[cipaordv1_fItop])*STATES[cipaordv1_a]*ALGEBRAIC[cipaordv1_i]+ ALGEBRAIC[cipaordv1_fItop]*STATES[cipaordv1_ap]*ALGEBRAIC[cipaordv1_ip]);
ALGEBRAIC[cipaordv1_IKr] =  CONSTANTS[cipaordv1_GKr]* pow((CONSTANTS[cipaordv1_ko]/5.40000), 1.0 / 2)*STATES[cipaordv1_O]*(STATES[cipaordv1_V] - ALGEBRAIC[cipaordv1_EK]);
ALGEBRAIC[cipaordv1_EKs] =  (( CONSTANTS[cipaordv1_R]*CONSTANTS[cipaordv1_T])/CONSTANTS[cipaordv1_F])*log((CONSTANTS[cipaordv1_ko]+ CONSTANTS[cipaordv1_PKNa]*CONSTANTS[cipaordv1_nao])/(STATES[cipaordv1_ki]+ CONSTANTS[cipaordv1_PKNa]*STATES[cipaordv1_nai]));
ALGEBRAIC[cipaordv1_KsCa] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[cipaordv1_cai], 1.40000));
ALGEBRAIC[cipaordv1_IKs] =  CONSTANTS[cipaordv1_GKs]*ALGEBRAIC[cipaordv1_KsCa]*STATES[cipaordv1_xs1]*STATES[cipaordv1_xs2]*(STATES[cipaordv1_V] - ALGEBRAIC[cipaordv1_EKs]);
ALGEBRAIC[cipaordv1_rk1] = 1.00000/(1.00000+exp(((STATES[cipaordv1_V]+105.800) -  2.60000*CONSTANTS[cipaordv1_ko])/9.49300));
ALGEBRAIC[cipaordv1_IK1] =  CONSTANTS[cipaordv1_GK1]* pow(CONSTANTS[cipaordv1_ko], 1.0 / 2)*ALGEBRAIC[cipaordv1_rk1]*STATES[cipaordv1_xk1]*(STATES[cipaordv1_V] - ALGEBRAIC[cipaordv1_EK]);
ALGEBRAIC[cipaordv1_Knao] =  CONSTANTS[cipaordv1_Knao0]*exp(( (1.00000 - CONSTANTS[cipaordv1_delta])*STATES[cipaordv1_V]*CONSTANTS[cipaordv1_F])/( 3.00000*CONSTANTS[cipaordv1_R]*CONSTANTS[cipaordv1_T]));
ALGEBRAIC[cipaordv1_a3] = ( CONSTANTS[cipaordv1_k3p]*pow(CONSTANTS[cipaordv1_ko]/CONSTANTS[cipaordv1_Kko], 2.00000))/((pow(1.00000+CONSTANTS[cipaordv1_nao]/ALGEBRAIC[cipaordv1_Knao], 3.00000)+pow(1.00000+CONSTANTS[cipaordv1_ko]/CONSTANTS[cipaordv1_Kko], 2.00000)) - 1.00000);
ALGEBRAIC[cipaordv1_P] = CONSTANTS[cipaordv1_eP]/(1.00000+CONSTANTS[cipaordv1_H]/CONSTANTS[cipaordv1_Khp]+STATES[cipaordv1_nai]/CONSTANTS[cipaordv1_Knap]+STATES[cipaordv1_ki]/CONSTANTS[cipaordv1_Kxkur]);
ALGEBRAIC[cipaordv1_b3] = ( CONSTANTS[cipaordv1_k3m]*ALGEBRAIC[cipaordv1_P]*CONSTANTS[cipaordv1_H])/(1.00000+CONSTANTS[cipaordv1_MgATP]/CONSTANTS[cipaordv1_Kmgatp]);
ALGEBRAIC[cipaordv1_Knai] =  CONSTANTS[cipaordv1_Knai0]*exp(( CONSTANTS[cipaordv1_delta]*STATES[cipaordv1_V]*CONSTANTS[cipaordv1_F])/( 3.00000*CONSTANTS[cipaordv1_R]*CONSTANTS[cipaordv1_T]));
ALGEBRAIC[cipaordv1_a1] = ( CONSTANTS[cipaordv1_k1p]*pow(STATES[cipaordv1_nai]/ALGEBRAIC[cipaordv1_Knai], 3.00000))/((pow(1.00000+STATES[cipaordv1_nai]/ALGEBRAIC[cipaordv1_Knai], 3.00000)+pow(1.00000+STATES[cipaordv1_ki]/CONSTANTS[cipaordv1_Kki], 2.00000)) - 1.00000);
ALGEBRAIC[cipaordv1_b2] = ( CONSTANTS[cipaordv1_k2m]*pow(CONSTANTS[cipaordv1_nao]/ALGEBRAIC[cipaordv1_Knao], 3.00000))/((pow(1.00000+CONSTANTS[cipaordv1_nao]/ALGEBRAIC[cipaordv1_Knao], 3.00000)+pow(1.00000+CONSTANTS[cipaordv1_ko]/CONSTANTS[cipaordv1_Kko], 2.00000)) - 1.00000);
ALGEBRAIC[cipaordv1_b4] = ( CONSTANTS[cipaordv1_k4m]*pow(STATES[cipaordv1_ki]/CONSTANTS[cipaordv1_Kki], 2.00000))/((pow(1.00000+STATES[cipaordv1_nai]/ALGEBRAIC[cipaordv1_Knai], 3.00000)+pow(1.00000+STATES[cipaordv1_ki]/CONSTANTS[cipaordv1_Kki], 2.00000)) - 1.00000);
ALGEBRAIC[cipaordv1_x1] =  CONSTANTS[cipaordv1_a4]*ALGEBRAIC[cipaordv1_a1]*CONSTANTS[cipaordv1_a2]+ ALGEBRAIC[cipaordv1_b2]*ALGEBRAIC[cipaordv1_b4]*ALGEBRAIC[cipaordv1_b3]+ CONSTANTS[cipaordv1_a2]*ALGEBRAIC[cipaordv1_b4]*ALGEBRAIC[cipaordv1_b3]+ ALGEBRAIC[cipaordv1_b3]*ALGEBRAIC[cipaordv1_a1]*CONSTANTS[cipaordv1_a2];
ALGEBRAIC[cipaordv1_x2] =  ALGEBRAIC[cipaordv1_b2]*CONSTANTS[cipaordv1_b1]*ALGEBRAIC[cipaordv1_b4]+ ALGEBRAIC[cipaordv1_a1]*CONSTANTS[cipaordv1_a2]*ALGEBRAIC[cipaordv1_a3]+ ALGEBRAIC[cipaordv1_a3]*CONSTANTS[cipaordv1_b1]*ALGEBRAIC[cipaordv1_b4]+ CONSTANTS[cipaordv1_a2]*ALGEBRAIC[cipaordv1_a3]*ALGEBRAIC[cipaordv1_b4];
ALGEBRAIC[cipaordv1_x3] =  CONSTANTS[cipaordv1_a2]*ALGEBRAIC[cipaordv1_a3]*CONSTANTS[cipaordv1_a4]+ ALGEBRAIC[cipaordv1_b3]*ALGEBRAIC[cipaordv1_b2]*CONSTANTS[cipaordv1_b1]+ ALGEBRAIC[cipaordv1_b2]*CONSTANTS[cipaordv1_b1]*CONSTANTS[cipaordv1_a4]+ ALGEBRAIC[cipaordv1_a3]*CONSTANTS[cipaordv1_a4]*CONSTANTS[cipaordv1_b1];
ALGEBRAIC[cipaordv1_x4] =  ALGEBRAIC[cipaordv1_b4]*ALGEBRAIC[cipaordv1_b3]*ALGEBRAIC[cipaordv1_b2]+ ALGEBRAIC[cipaordv1_a3]*CONSTANTS[cipaordv1_a4]*ALGEBRAIC[cipaordv1_a1]+ ALGEBRAIC[cipaordv1_b2]*CONSTANTS[cipaordv1_a4]*ALGEBRAIC[cipaordv1_a1]+ ALGEBRAIC[cipaordv1_b3]*ALGEBRAIC[cipaordv1_b2]*ALGEBRAIC[cipaordv1_a1];
ALGEBRAIC[cipaordv1_E1] = ALGEBRAIC[cipaordv1_x1]/(ALGEBRAIC[cipaordv1_x1]+ALGEBRAIC[cipaordv1_x2]+ALGEBRAIC[cipaordv1_x3]+ALGEBRAIC[cipaordv1_x4]);
ALGEBRAIC[cipaordv1_E2] = ALGEBRAIC[cipaordv1_x2]/(ALGEBRAIC[cipaordv1_x1]+ALGEBRAIC[cipaordv1_x2]+ALGEBRAIC[cipaordv1_x3]+ALGEBRAIC[cipaordv1_x4]);
ALGEBRAIC[cipaordv1_JnakNa] =  3.00000*( ALGEBRAIC[cipaordv1_E1]*ALGEBRAIC[cipaordv1_a3] -  ALGEBRAIC[cipaordv1_E2]*ALGEBRAIC[cipaordv1_b3]);
ALGEBRAIC[cipaordv1_E3] = ALGEBRAIC[cipaordv1_x3]/(ALGEBRAIC[cipaordv1_x1]+ALGEBRAIC[cipaordv1_x2]+ALGEBRAIC[cipaordv1_x3]+ALGEBRAIC[cipaordv1_x4]);
ALGEBRAIC[cipaordv1_E4] = ALGEBRAIC[cipaordv1_x4]/(ALGEBRAIC[cipaordv1_x1]+ALGEBRAIC[cipaordv1_x2]+ALGEBRAIC[cipaordv1_x3]+ALGEBRAIC[cipaordv1_x4]);
ALGEBRAIC[cipaordv1_JnakK] =  2.00000*( ALGEBRAIC[cipaordv1_E4]*CONSTANTS[cipaordv1_b1] -  ALGEBRAIC[cipaordv1_E3]*ALGEBRAIC[cipaordv1_a1]);
ALGEBRAIC[cipaordv1_INaK] =  CONSTANTS[cipaordv1_Pnak]*( CONSTANTS[cipaordv1_zna]*ALGEBRAIC[cipaordv1_JnakNa]+ CONSTANTS[cipaordv1_zk]*ALGEBRAIC[cipaordv1_JnakK]);
ALGEBRAIC[cipaordv1_xkb] = 1.00000/(1.00000+exp(- (STATES[cipaordv1_V] - 14.4800)/18.3400));
ALGEBRAIC[cipaordv1_IKb] =  CONSTANTS[cipaordv1_GKb]*ALGEBRAIC[cipaordv1_xkb]*(STATES[cipaordv1_V] - ALGEBRAIC[cipaordv1_EK]);
ALGEBRAIC[cipaordv1_Istim] = (TIME>=CONSTANTS[cipaordv1_stim_start]&&TIME<=CONSTANTS[cipaordv1_stim_end]&&(TIME - CONSTANTS[cipaordv1_stim_start]) -  floor((TIME - CONSTANTS[cipaordv1_stim_start])/CONSTANTS[cipaordv1_BCL])*CONSTANTS[cipaordv1_BCL]<=CONSTANTS[cipaordv1_duration] ? CONSTANTS[cipaordv1_amp] : 0.000000);
ALGEBRAIC[cipaordv1_JdiffK] = (STATES[cipaordv1_kss] - STATES[cipaordv1_ki])/2.00000;
ALGEBRAIC[cipaordv1_A_3] = ( 0.750000*CONSTANTS[cipaordv1_ffrt]*( STATES[cipaordv1_kss]*exp(ALGEBRAIC[cipaordv1_vfrt]) - CONSTANTS[cipaordv1_ko]))/CONSTANTS[cipaordv1_B_3];
ALGEBRAIC[cipaordv1_U_3] =  CONSTANTS[cipaordv1_B_3]*(STATES[cipaordv1_V] - CONSTANTS[cipaordv1_v0_CaL]);
ALGEBRAIC[cipaordv1_PhiCaK] = (- 1.00000e-07<=ALGEBRAIC[cipaordv1_U_3]&&ALGEBRAIC[cipaordv1_U_3]<=1.00000e-07 ?  ALGEBRAIC[cipaordv1_A_3]*(1.00000 -  0.500000*ALGEBRAIC[cipaordv1_U_3]) : ( ALGEBRAIC[cipaordv1_A_3]*ALGEBRAIC[cipaordv1_U_3])/(exp(ALGEBRAIC[cipaordv1_U_3]) - 1.00000));
ALGEBRAIC[cipaordv1_ICaK] =  (1.00000 - ALGEBRAIC[cipaordv1_fICaLp])*CONSTANTS[cipaordv1_PCaK]*ALGEBRAIC[cipaordv1_PhiCaK]*STATES[cipaordv1_d]*( ALGEBRAIC[cipaordv1_f]*(1.00000 - STATES[cipaordv1_nca])+ STATES[cipaordv1_jca]*ALGEBRAIC[cipaordv1_fca]*STATES[cipaordv1_nca])+ ALGEBRAIC[cipaordv1_fICaLp]*CONSTANTS[cipaordv1_PCaKp]*ALGEBRAIC[cipaordv1_PhiCaK]*STATES[cipaordv1_d]*( ALGEBRAIC[cipaordv1_fp]*(1.00000 - STATES[cipaordv1_nca])+ STATES[cipaordv1_jca]*ALGEBRAIC[cipaordv1_fcap]*STATES[cipaordv1_nca]);
ALGEBRAIC[cipaordv1_ENa] =  (( CONSTANTS[cipaordv1_R]*CONSTANTS[cipaordv1_T])/CONSTANTS[cipaordv1_F])*log(CONSTANTS[cipaordv1_nao]/STATES[cipaordv1_nai]);
ALGEBRAIC[cipaordv1_h] =  CONSTANTS[cipaordv1_Ahf]*STATES[cipaordv1_hf]+ CONSTANTS[cipaordv1_Ahs]*STATES[cipaordv1_hs];
ALGEBRAIC[cipaordv1_hp] =  CONSTANTS[cipaordv1_Ahf]*STATES[cipaordv1_hf]+ CONSTANTS[cipaordv1_Ahs]*STATES[cipaordv1_hsp];
ALGEBRAIC[cipaordv1_fINap] = 1.00000/(1.00000+CONSTANTS[cipaordv1_KmCaMK]/ALGEBRAIC[cipaordv1_CaMKa]);
ALGEBRAIC[cipaordv1_INa] =  CONSTANTS[cipaordv1_GNa]*(STATES[cipaordv1_V] - ALGEBRAIC[cipaordv1_ENa])*pow(STATES[cipaordv1_m], 3.00000)*( (1.00000 - ALGEBRAIC[cipaordv1_fINap])*ALGEBRAIC[cipaordv1_h]*STATES[cipaordv1_j]+ ALGEBRAIC[cipaordv1_fINap]*ALGEBRAIC[cipaordv1_hp]*STATES[cipaordv1_jp]);
ALGEBRAIC[cipaordv1_fINaLp] = 1.00000/(1.00000+CONSTANTS[cipaordv1_KmCaMK]/ALGEBRAIC[cipaordv1_CaMKa]);
ALGEBRAIC[cipaordv1_INaL] =  CONSTANTS[cipaordv1_GNaL]*(STATES[cipaordv1_V] - ALGEBRAIC[cipaordv1_ENa])*STATES[cipaordv1_mL]*( (1.00000 - ALGEBRAIC[cipaordv1_fINaLp])*STATES[cipaordv1_hL]+ ALGEBRAIC[cipaordv1_fINaLp]*STATES[cipaordv1_hLp]);
ALGEBRAIC[cipaordv1_allo_i] = 1.00000/(1.00000+pow(CONSTANTS[cipaordv1_KmCaAct]/STATES[cipaordv1_cai], 2.00000));
ALGEBRAIC[cipaordv1_hna] = exp(( CONSTANTS[cipaordv1_qna]*STATES[cipaordv1_V]*CONSTANTS[cipaordv1_F])/( CONSTANTS[cipaordv1_R]*CONSTANTS[cipaordv1_T]));
ALGEBRAIC[cipaordv1_h7_i] = 1.00000+ (CONSTANTS[cipaordv1_nao]/CONSTANTS[cipaordv1_kna3])*(1.00000+1.00000/ALGEBRAIC[cipaordv1_hna]);
ALGEBRAIC[cipaordv1_h8_i] = CONSTANTS[cipaordv1_nao]/( CONSTANTS[cipaordv1_kna3]*ALGEBRAIC[cipaordv1_hna]*ALGEBRAIC[cipaordv1_h7_i]);
ALGEBRAIC[cipaordv1_k3pp_i] =  ALGEBRAIC[cipaordv1_h8_i]*CONSTANTS[cipaordv1_wnaca];
ALGEBRAIC[cipaordv1_h1_i] = 1.00000+ (STATES[cipaordv1_nai]/CONSTANTS[cipaordv1_kna3])*(1.00000+ALGEBRAIC[cipaordv1_hna]);
ALGEBRAIC[cipaordv1_h2_i] = ( STATES[cipaordv1_nai]*ALGEBRAIC[cipaordv1_hna])/( CONSTANTS[cipaordv1_kna3]*ALGEBRAIC[cipaordv1_h1_i]);
ALGEBRAIC[cipaordv1_k4pp_i] =  ALGEBRAIC[cipaordv1_h2_i]*CONSTANTS[cipaordv1_wnaca];
ALGEBRAIC[cipaordv1_h4_i] = 1.00000+ (STATES[cipaordv1_nai]/CONSTANTS[cipaordv1_kna1])*(1.00000+STATES[cipaordv1_nai]/CONSTANTS[cipaordv1_kna2]);
ALGEBRAIC[cipaordv1_h5_i] = ( STATES[cipaordv1_nai]*STATES[cipaordv1_nai])/( ALGEBRAIC[cipaordv1_h4_i]*CONSTANTS[cipaordv1_kna1]*CONSTANTS[cipaordv1_kna2]);
ALGEBRAIC[cipaordv1_k7_i] =  ALGEBRAIC[cipaordv1_h5_i]*ALGEBRAIC[cipaordv1_h2_i]*CONSTANTS[cipaordv1_wna];
ALGEBRAIC[cipaordv1_k8_i] =  ALGEBRAIC[cipaordv1_h8_i]*CONSTANTS[cipaordv1_h11_i]*CONSTANTS[cipaordv1_wna];
ALGEBRAIC[cipaordv1_h9_i] = 1.00000/ALGEBRAIC[cipaordv1_h7_i];
ALGEBRAIC[cipaordv1_k3p_i] =  ALGEBRAIC[cipaordv1_h9_i]*CONSTANTS[cipaordv1_wca];
ALGEBRAIC[cipaordv1_k3_i] = ALGEBRAIC[cipaordv1_k3p_i]+ALGEBRAIC[cipaordv1_k3pp_i];
ALGEBRAIC[cipaordv1_hca] = exp(( CONSTANTS[cipaordv1_qca]*STATES[cipaordv1_V]*CONSTANTS[cipaordv1_F])/( CONSTANTS[cipaordv1_R]*CONSTANTS[cipaordv1_T]));
ALGEBRAIC[cipaordv1_h3_i] = 1.00000/ALGEBRAIC[cipaordv1_h1_i];
ALGEBRAIC[cipaordv1_k4p_i] = ( ALGEBRAIC[cipaordv1_h3_i]*CONSTANTS[cipaordv1_wca])/ALGEBRAIC[cipaordv1_hca];
ALGEBRAIC[cipaordv1_k4_i] = ALGEBRAIC[cipaordv1_k4p_i]+ALGEBRAIC[cipaordv1_k4pp_i];
ALGEBRAIC[cipaordv1_h6_i] = 1.00000/ALGEBRAIC[cipaordv1_h4_i];
ALGEBRAIC[cipaordv1_k6_i] =  ALGEBRAIC[cipaordv1_h6_i]*STATES[cipaordv1_cai]*CONSTANTS[cipaordv1_kcaon];
ALGEBRAIC[cipaordv1_x1_i] =  CONSTANTS[cipaordv1_k2_i]*ALGEBRAIC[cipaordv1_k4_i]*(ALGEBRAIC[cipaordv1_k7_i]+ALGEBRAIC[cipaordv1_k6_i])+ CONSTANTS[cipaordv1_k5_i]*ALGEBRAIC[cipaordv1_k7_i]*(CONSTANTS[cipaordv1_k2_i]+ALGEBRAIC[cipaordv1_k3_i]);
ALGEBRAIC[cipaordv1_x2_i] =  CONSTANTS[cipaordv1_k1_i]*ALGEBRAIC[cipaordv1_k7_i]*(ALGEBRAIC[cipaordv1_k4_i]+CONSTANTS[cipaordv1_k5_i])+ ALGEBRAIC[cipaordv1_k4_i]*ALGEBRAIC[cipaordv1_k6_i]*(CONSTANTS[cipaordv1_k1_i]+ALGEBRAIC[cipaordv1_k8_i]);
ALGEBRAIC[cipaordv1_x3_i] =  CONSTANTS[cipaordv1_k1_i]*ALGEBRAIC[cipaordv1_k3_i]*(ALGEBRAIC[cipaordv1_k7_i]+ALGEBRAIC[cipaordv1_k6_i])+ ALGEBRAIC[cipaordv1_k8_i]*ALGEBRAIC[cipaordv1_k6_i]*(CONSTANTS[cipaordv1_k2_i]+ALGEBRAIC[cipaordv1_k3_i]);
ALGEBRAIC[cipaordv1_x4_i] =  CONSTANTS[cipaordv1_k2_i]*ALGEBRAIC[cipaordv1_k8_i]*(ALGEBRAIC[cipaordv1_k4_i]+CONSTANTS[cipaordv1_k5_i])+ ALGEBRAIC[cipaordv1_k3_i]*CONSTANTS[cipaordv1_k5_i]*(CONSTANTS[cipaordv1_k1_i]+ALGEBRAIC[cipaordv1_k8_i]);
ALGEBRAIC[cipaordv1_E1_i] = ALGEBRAIC[cipaordv1_x1_i]/(ALGEBRAIC[cipaordv1_x1_i]+ALGEBRAIC[cipaordv1_x2_i]+ALGEBRAIC[cipaordv1_x3_i]+ALGEBRAIC[cipaordv1_x4_i]);
ALGEBRAIC[cipaordv1_E2_i] = ALGEBRAIC[cipaordv1_x2_i]/(ALGEBRAIC[cipaordv1_x1_i]+ALGEBRAIC[cipaordv1_x2_i]+ALGEBRAIC[cipaordv1_x3_i]+ALGEBRAIC[cipaordv1_x4_i]);
ALGEBRAIC[cipaordv1_E3_i] = ALGEBRAIC[cipaordv1_x3_i]/(ALGEBRAIC[cipaordv1_x1_i]+ALGEBRAIC[cipaordv1_x2_i]+ALGEBRAIC[cipaordv1_x3_i]+ALGEBRAIC[cipaordv1_x4_i]);
ALGEBRAIC[cipaordv1_E4_i] = ALGEBRAIC[cipaordv1_x4_i]/(ALGEBRAIC[cipaordv1_x1_i]+ALGEBRAIC[cipaordv1_x2_i]+ALGEBRAIC[cipaordv1_x3_i]+ALGEBRAIC[cipaordv1_x4_i]);
ALGEBRAIC[cipaordv1_JncxNa_i] = ( 3.00000*( ALGEBRAIC[cipaordv1_E4_i]*ALGEBRAIC[cipaordv1_k7_i] -  ALGEBRAIC[cipaordv1_E1_i]*ALGEBRAIC[cipaordv1_k8_i])+ ALGEBRAIC[cipaordv1_E3_i]*ALGEBRAIC[cipaordv1_k4pp_i]) -  ALGEBRAIC[cipaordv1_E2_i]*ALGEBRAIC[cipaordv1_k3pp_i];
ALGEBRAIC[cipaordv1_JncxCa_i] =  ALGEBRAIC[cipaordv1_E2_i]*CONSTANTS[cipaordv1_k2_i] -  ALGEBRAIC[cipaordv1_E1_i]*CONSTANTS[cipaordv1_k1_i];
ALGEBRAIC[cipaordv1_INaCa_i] =  0.800000*CONSTANTS[cipaordv1_Gncx]*ALGEBRAIC[cipaordv1_allo_i]*( CONSTANTS[cipaordv1_zna]*ALGEBRAIC[cipaordv1_JncxNa_i]+ CONSTANTS[cipaordv1_zca]*ALGEBRAIC[cipaordv1_JncxCa_i]);
ALGEBRAIC[cipaordv1_A_Nab] = ( CONSTANTS[cipaordv1_PNab]*CONSTANTS[cipaordv1_ffrt]*( STATES[cipaordv1_nai]*exp(ALGEBRAIC[cipaordv1_vfrt]) - CONSTANTS[cipaordv1_nao]))/CONSTANTS[cipaordv1_B_Nab];
ALGEBRAIC[cipaordv1_U_Nab] =  CONSTANTS[cipaordv1_B_Nab]*(STATES[cipaordv1_V] - CONSTANTS[cipaordv1_v0_Nab]);
ALGEBRAIC[cipaordv1_INab] = (- 1.00000e-07<=ALGEBRAIC[cipaordv1_U_Nab]&&ALGEBRAIC[cipaordv1_U_Nab]<=1.00000e-07 ?  ALGEBRAIC[cipaordv1_A_Nab]*(1.00000 -  0.500000*ALGEBRAIC[cipaordv1_U_Nab]) : ( ALGEBRAIC[cipaordv1_A_Nab]*ALGEBRAIC[cipaordv1_U_Nab])/(exp(ALGEBRAIC[cipaordv1_U_Nab]) - 1.00000));
ALGEBRAIC[cipaordv1_JdiffNa] = (STATES[cipaordv1_nass] - STATES[cipaordv1_nai])/2.00000;
ALGEBRAIC[cipaordv1_A_2] = ( 0.750000*CONSTANTS[cipaordv1_ffrt]*( STATES[cipaordv1_nass]*exp(ALGEBRAIC[cipaordv1_vfrt]) - CONSTANTS[cipaordv1_nao]))/CONSTANTS[cipaordv1_B_2];
ALGEBRAIC[cipaordv1_U_2] =  CONSTANTS[cipaordv1_B_2]*(STATES[cipaordv1_V] - CONSTANTS[cipaordv1_v0_CaL]);
ALGEBRAIC[cipaordv1_PhiCaNa] = (- 1.00000e-07<=ALGEBRAIC[cipaordv1_U_2]&&ALGEBRAIC[cipaordv1_U_2]<=1.00000e-07 ?  ALGEBRAIC[cipaordv1_A_2]*(1.00000 -  0.500000*ALGEBRAIC[cipaordv1_U_2]) : ( ALGEBRAIC[cipaordv1_A_2]*ALGEBRAIC[cipaordv1_U_2])/(exp(ALGEBRAIC[cipaordv1_U_2]) - 1.00000));
ALGEBRAIC[cipaordv1_ICaNa] =  (1.00000 - ALGEBRAIC[cipaordv1_fICaLp])*CONSTANTS[cipaordv1_PCaNa]*ALGEBRAIC[cipaordv1_PhiCaNa]*STATES[cipaordv1_d]*( ALGEBRAIC[cipaordv1_f]*(1.00000 - STATES[cipaordv1_nca])+ STATES[cipaordv1_jca]*ALGEBRAIC[cipaordv1_fca]*STATES[cipaordv1_nca])+ ALGEBRAIC[cipaordv1_fICaLp]*CONSTANTS[cipaordv1_PCaNap]*ALGEBRAIC[cipaordv1_PhiCaNa]*STATES[cipaordv1_d]*( ALGEBRAIC[cipaordv1_fp]*(1.00000 - STATES[cipaordv1_nca])+ STATES[cipaordv1_jca]*ALGEBRAIC[cipaordv1_fcap]*STATES[cipaordv1_nca]);
ALGEBRAIC[cipaordv1_allo_ss] = 1.00000/(1.00000+pow(CONSTANTS[cipaordv1_KmCaAct]/STATES[cipaordv1_cass], 2.00000));
ALGEBRAIC[cipaordv1_h7_ss] = 1.00000+ (CONSTANTS[cipaordv1_nao]/CONSTANTS[cipaordv1_kna3])*(1.00000+1.00000/ALGEBRAIC[cipaordv1_hna]);
ALGEBRAIC[cipaordv1_h8_ss] = CONSTANTS[cipaordv1_nao]/( CONSTANTS[cipaordv1_kna3]*ALGEBRAIC[cipaordv1_hna]*ALGEBRAIC[cipaordv1_h7_ss]);
ALGEBRAIC[cipaordv1_k3pp_ss] =  ALGEBRAIC[cipaordv1_h8_ss]*CONSTANTS[cipaordv1_wnaca];
ALGEBRAIC[cipaordv1_h1_ss] = 1.00000+ (STATES[cipaordv1_nass]/CONSTANTS[cipaordv1_kna3])*(1.00000+ALGEBRAIC[cipaordv1_hna]);
ALGEBRAIC[cipaordv1_h2_ss] = ( STATES[cipaordv1_nass]*ALGEBRAIC[cipaordv1_hna])/( CONSTANTS[cipaordv1_kna3]*ALGEBRAIC[cipaordv1_h1_ss]);
ALGEBRAIC[cipaordv1_k4pp_ss] =  ALGEBRAIC[cipaordv1_h2_ss]*CONSTANTS[cipaordv1_wnaca];
ALGEBRAIC[cipaordv1_h4_ss] = 1.00000+ (STATES[cipaordv1_nass]/CONSTANTS[cipaordv1_kna1])*(1.00000+STATES[cipaordv1_nass]/CONSTANTS[cipaordv1_kna2]);
ALGEBRAIC[cipaordv1_h5_ss] = ( STATES[cipaordv1_nass]*STATES[cipaordv1_nass])/( ALGEBRAIC[cipaordv1_h4_ss]*CONSTANTS[cipaordv1_kna1]*CONSTANTS[cipaordv1_kna2]);
ALGEBRAIC[cipaordv1_k7_ss] =  ALGEBRAIC[cipaordv1_h5_ss]*ALGEBRAIC[cipaordv1_h2_ss]*CONSTANTS[cipaordv1_wna];
ALGEBRAIC[cipaordv1_k8_ss] =  ALGEBRAIC[cipaordv1_h8_ss]*CONSTANTS[cipaordv1_h11_ss]*CONSTANTS[cipaordv1_wna];
ALGEBRAIC[cipaordv1_h9_ss] = 1.00000/ALGEBRAIC[cipaordv1_h7_ss];
ALGEBRAIC[cipaordv1_k3p_ss] =  ALGEBRAIC[cipaordv1_h9_ss]*CONSTANTS[cipaordv1_wca];
ALGEBRAIC[cipaordv1_k3_ss] = ALGEBRAIC[cipaordv1_k3p_ss]+ALGEBRAIC[cipaordv1_k3pp_ss];
ALGEBRAIC[cipaordv1_h3_ss] = 1.00000/ALGEBRAIC[cipaordv1_h1_ss];
ALGEBRAIC[cipaordv1_k4p_ss] = ( ALGEBRAIC[cipaordv1_h3_ss]*CONSTANTS[cipaordv1_wca])/ALGEBRAIC[cipaordv1_hca];
ALGEBRAIC[cipaordv1_k4_ss] = ALGEBRAIC[cipaordv1_k4p_ss]+ALGEBRAIC[cipaordv1_k4pp_ss];
ALGEBRAIC[cipaordv1_h6_ss] = 1.00000/ALGEBRAIC[cipaordv1_h4_ss];
ALGEBRAIC[cipaordv1_k6_ss] =  ALGEBRAIC[cipaordv1_h6_ss]*STATES[cipaordv1_cass]*CONSTANTS[cipaordv1_kcaon];
ALGEBRAIC[cipaordv1_x1_ss] =  CONSTANTS[cipaordv1_k2_ss]*ALGEBRAIC[cipaordv1_k4_ss]*(ALGEBRAIC[cipaordv1_k7_ss]+ALGEBRAIC[cipaordv1_k6_ss])+ CONSTANTS[cipaordv1_k5_ss]*ALGEBRAIC[cipaordv1_k7_ss]*(CONSTANTS[cipaordv1_k2_ss]+ALGEBRAIC[cipaordv1_k3_ss]);
ALGEBRAIC[cipaordv1_x2_ss] =  CONSTANTS[cipaordv1_k1_ss]*ALGEBRAIC[cipaordv1_k7_ss]*(ALGEBRAIC[cipaordv1_k4_ss]+CONSTANTS[cipaordv1_k5_ss])+ ALGEBRAIC[cipaordv1_k4_ss]*ALGEBRAIC[cipaordv1_k6_ss]*(CONSTANTS[cipaordv1_k1_ss]+ALGEBRAIC[cipaordv1_k8_ss]);
ALGEBRAIC[cipaordv1_x3_ss] =  CONSTANTS[cipaordv1_k1_ss]*ALGEBRAIC[cipaordv1_k3_ss]*(ALGEBRAIC[cipaordv1_k7_ss]+ALGEBRAIC[cipaordv1_k6_ss])+ ALGEBRAIC[cipaordv1_k8_ss]*ALGEBRAIC[cipaordv1_k6_ss]*(CONSTANTS[cipaordv1_k2_ss]+ALGEBRAIC[cipaordv1_k3_ss]);
ALGEBRAIC[cipaordv1_x4_ss] =  CONSTANTS[cipaordv1_k2_ss]*ALGEBRAIC[cipaordv1_k8_ss]*(ALGEBRAIC[cipaordv1_k4_ss]+CONSTANTS[cipaordv1_k5_ss])+ ALGEBRAIC[cipaordv1_k3_ss]*CONSTANTS[cipaordv1_k5_ss]*(CONSTANTS[cipaordv1_k1_ss]+ALGEBRAIC[cipaordv1_k8_ss]);
ALGEBRAIC[cipaordv1_E1_ss] = ALGEBRAIC[cipaordv1_x1_ss]/(ALGEBRAIC[cipaordv1_x1_ss]+ALGEBRAIC[cipaordv1_x2_ss]+ALGEBRAIC[cipaordv1_x3_ss]+ALGEBRAIC[cipaordv1_x4_ss]);
ALGEBRAIC[cipaordv1_E2_ss] = ALGEBRAIC[cipaordv1_x2_ss]/(ALGEBRAIC[cipaordv1_x1_ss]+ALGEBRAIC[cipaordv1_x2_ss]+ALGEBRAIC[cipaordv1_x3_ss]+ALGEBRAIC[cipaordv1_x4_ss]);
ALGEBRAIC[cipaordv1_E3_ss] = ALGEBRAIC[cipaordv1_x3_ss]/(ALGEBRAIC[cipaordv1_x1_ss]+ALGEBRAIC[cipaordv1_x2_ss]+ALGEBRAIC[cipaordv1_x3_ss]+ALGEBRAIC[cipaordv1_x4_ss]);
ALGEBRAIC[cipaordv1_E4_ss] = ALGEBRAIC[cipaordv1_x4_ss]/(ALGEBRAIC[cipaordv1_x1_ss]+ALGEBRAIC[cipaordv1_x2_ss]+ALGEBRAIC[cipaordv1_x3_ss]+ALGEBRAIC[cipaordv1_x4_ss]);
ALGEBRAIC[cipaordv1_JncxNa_ss] = ( 3.00000*( ALGEBRAIC[cipaordv1_E4_ss]*ALGEBRAIC[cipaordv1_k7_ss] -  ALGEBRAIC[cipaordv1_E1_ss]*ALGEBRAIC[cipaordv1_k8_ss])+ ALGEBRAIC[cipaordv1_E3_ss]*ALGEBRAIC[cipaordv1_k4pp_ss]) -  ALGEBRAIC[cipaordv1_E2_ss]*ALGEBRAIC[cipaordv1_k3pp_ss];
ALGEBRAIC[cipaordv1_JncxCa_ss] =  ALGEBRAIC[cipaordv1_E2_ss]*CONSTANTS[cipaordv1_k2_ss] -  ALGEBRAIC[cipaordv1_E1_ss]*CONSTANTS[cipaordv1_k1_ss];
ALGEBRAIC[cipaordv1_INaCa_ss] =  0.200000*CONSTANTS[cipaordv1_Gncx]*ALGEBRAIC[cipaordv1_allo_ss]*( CONSTANTS[cipaordv1_zna]*ALGEBRAIC[cipaordv1_JncxNa_ss]+ CONSTANTS[cipaordv1_zca]*ALGEBRAIC[cipaordv1_JncxCa_ss]);
ALGEBRAIC[cipaordv1_IpCa] = ( CONSTANTS[cipaordv1_GpCa]*STATES[cipaordv1_cai])/(CONSTANTS[cipaordv1_KmCap]+STATES[cipaordv1_cai]);
ALGEBRAIC[cipaordv1_A_Cab] = ( CONSTANTS[cipaordv1_PCab]*4.00000*CONSTANTS[cipaordv1_ffrt]*( STATES[cipaordv1_cai]*exp( 2.00000*ALGEBRAIC[cipaordv1_vfrt]) -  0.341000*CONSTANTS[cipaordv1_cao]))/CONSTANTS[cipaordv1_B_Cab];
ALGEBRAIC[cipaordv1_U_Cab] =  CONSTANTS[cipaordv1_B_Cab]*(STATES[cipaordv1_V] - CONSTANTS[cipaordv1_v0_Cab]);
ALGEBRAIC[cipaordv1_ICab] = (- 1.00000e-07<=ALGEBRAIC[cipaordv1_U_Cab]&&ALGEBRAIC[cipaordv1_U_Cab]<=1.00000e-07 ?  ALGEBRAIC[cipaordv1_A_Cab]*(1.00000 -  0.500000*ALGEBRAIC[cipaordv1_U_Cab]) : ( ALGEBRAIC[cipaordv1_A_Cab]*ALGEBRAIC[cipaordv1_U_Cab])/(exp(ALGEBRAIC[cipaordv1_U_Cab]) - 1.00000));
ALGEBRAIC[cipaordv1_Jdiff] = (STATES[cipaordv1_cass] - STATES[cipaordv1_cai])/0.200000;
ALGEBRAIC[cipaordv1_fJrelp] = 1.00000/(1.00000+CONSTANTS[cipaordv1_KmCaMK]/ALGEBRAIC[cipaordv1_CaMKa]);
ALGEBRAIC[cipaordv1_Jrel] =  CONSTANTS[cipaordv1_Jrel_scaling_factor]*( (1.00000 - ALGEBRAIC[cipaordv1_fJrelp])*STATES[cipaordv1_Jrelnp]+ ALGEBRAIC[cipaordv1_fJrelp]*STATES[cipaordv1_Jrelp]);
ALGEBRAIC[cipaordv1_Bcass] = 1.00000/(1.00000+( CONSTANTS[cipaordv1_BSRmax]*CONSTANTS[cipaordv1_KmBSR])/pow(CONSTANTS[cipaordv1_KmBSR]+STATES[cipaordv1_cass], 2.00000)+( CONSTANTS[cipaordv1_BSLmax]*CONSTANTS[cipaordv1_KmBSL])/pow(CONSTANTS[cipaordv1_KmBSL]+STATES[cipaordv1_cass], 2.00000));
ALGEBRAIC[cipaordv1_Jupnp] = ( CONSTANTS[cipaordv1_upScale]*0.00437500*STATES[cipaordv1_cai])/(STATES[cipaordv1_cai]+0.000920000);
ALGEBRAIC[cipaordv1_Jupp] = ( CONSTANTS[cipaordv1_upScale]*2.75000*0.00437500*STATES[cipaordv1_cai])/((STATES[cipaordv1_cai]+0.000920000) - 0.000170000);
ALGEBRAIC[cipaordv1_fJupp] = 1.00000/(1.00000+CONSTANTS[cipaordv1_KmCaMK]/ALGEBRAIC[cipaordv1_CaMKa]);
ALGEBRAIC[cipaordv1_Jleak] = ( 0.00393750*STATES[cipaordv1_cansr])/15.0000;
ALGEBRAIC[cipaordv1_Jup] =  CONSTANTS[cipaordv1_Jup_b]*(( (1.00000 - ALGEBRAIC[cipaordv1_fJupp])*ALGEBRAIC[cipaordv1_Jupnp]+ ALGEBRAIC[cipaordv1_fJupp]*ALGEBRAIC[cipaordv1_Jupp]) - ALGEBRAIC[cipaordv1_Jleak]);
// ALGEBRAIC[cipaordv1_Bcai] = 1.00000/(1.00000+( CONSTANTS[cipaordv1_cmdnmax]*CONSTANTS[cipaordv1_kmcmdn])/pow(CONSTANTS[cipaordv1_kmcmdn]+STATES[cipaordv1_cai], 2.00000)+( CONSTANTS[cipaordv1_trpnmax]*CONSTANTS[cipaordv1_kmtrpn])/pow(CONSTANTS[cipaordv1_kmtrpn]+STATES[cipaordv1_cai], 2.00000));
ALGEBRAIC[cipaordv1_Bcai] = 1.00000/(1.00000+( CONSTANTS[cipaordv1_cmdnmax]*CONSTANTS[cipaordv1_kmcmdn])/pow(CONSTANTS[cipaordv1_kmcmdn]+STATES[cipaordv1_cai], 2.00000)); //modified
// using 2011's Bcai
ALGEBRAIC[cipaordv1_Jtr] = (STATES[cipaordv1_cansr] - STATES[cipaordv1_cajsr])/100.000;
ALGEBRAIC[cipaordv1_Bcajsr] = 1.00000/(1.00000+( CONSTANTS[cipaordv1_csqnmax]*CONSTANTS[cipaordv1_kmcsqn])/pow(CONSTANTS[cipaordv1_kmcsqn]+STATES[cipaordv1_cajsr], 2.00000));

//RATES[cipaordv1_D] = CONSTANTS[cipaordv1_cnc];
RATES[cipaordv1_D] = 0.;
RATES[cipaordv1_IC1] = (- ( CONSTANTS[cipaordv1_A11]*exp( CONSTANTS[cipaordv1_B11]*STATES[cipaordv1_V])*STATES[cipaordv1_IC1]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q11]))/10.0000) -  CONSTANTS[cipaordv1_A21]*exp( CONSTANTS[cipaordv1_B21]*STATES[cipaordv1_V])*STATES[cipaordv1_IC2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q21]))/10.0000))+ CONSTANTS[cipaordv1_A51]*exp( CONSTANTS[cipaordv1_B51]*STATES[cipaordv1_V])*STATES[cipaordv1_C1]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q51]))/10.0000)) -  CONSTANTS[cipaordv1_A61]*exp( CONSTANTS[cipaordv1_B61]*STATES[cipaordv1_V])*STATES[cipaordv1_IC1]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q61]))/10.0000);
RATES[cipaordv1_IC2] = ((( CONSTANTS[cipaordv1_A11]*exp( CONSTANTS[cipaordv1_B11]*STATES[cipaordv1_V])*STATES[cipaordv1_IC1]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q11]))/10.0000) -  CONSTANTS[cipaordv1_A21]*exp( CONSTANTS[cipaordv1_B21]*STATES[cipaordv1_V])*STATES[cipaordv1_IC2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q21]))/10.0000)) - ( CONSTANTS[cipaordv1_A3]*exp( CONSTANTS[cipaordv1_B3]*STATES[cipaordv1_V])*STATES[cipaordv1_IC2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q3]))/10.0000) -  CONSTANTS[cipaordv1_A4]*exp( CONSTANTS[cipaordv1_B4]*STATES[cipaordv1_V])*STATES[cipaordv1_IO]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q4]))/10.0000)))+ CONSTANTS[cipaordv1_A52]*exp( CONSTANTS[cipaordv1_B52]*STATES[cipaordv1_V])*STATES[cipaordv1_C2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q52]))/10.0000)) -  CONSTANTS[cipaordv1_A62]*exp( CONSTANTS[cipaordv1_B62]*STATES[cipaordv1_V])*STATES[cipaordv1_IC2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q62]))/10.0000);
RATES[cipaordv1_C1] = - ( CONSTANTS[cipaordv1_A1]*exp( CONSTANTS[cipaordv1_B1]*STATES[cipaordv1_V])*STATES[cipaordv1_C1]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q1]))/10.0000) -  CONSTANTS[cipaordv1_A2]*exp( CONSTANTS[cipaordv1_B2]*STATES[cipaordv1_V])*STATES[cipaordv1_C2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q2]))/10.0000)) - ( CONSTANTS[cipaordv1_A51]*exp( CONSTANTS[cipaordv1_B51]*STATES[cipaordv1_V])*STATES[cipaordv1_C1]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q51]))/10.0000) -  CONSTANTS[cipaordv1_A61]*exp( CONSTANTS[cipaordv1_B61]*STATES[cipaordv1_V])*STATES[cipaordv1_IC1]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q61]))/10.0000));
RATES[cipaordv1_C2] = (( CONSTANTS[cipaordv1_A1]*exp( CONSTANTS[cipaordv1_B1]*STATES[cipaordv1_V])*STATES[cipaordv1_C1]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q1]))/10.0000) -  CONSTANTS[cipaordv1_A2]*exp( CONSTANTS[cipaordv1_B2]*STATES[cipaordv1_V])*STATES[cipaordv1_C2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q2]))/10.0000)) - ( CONSTANTS[cipaordv1_A31]*exp( CONSTANTS[cipaordv1_B31]*STATES[cipaordv1_V])*STATES[cipaordv1_C2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q31]))/10.0000) -  CONSTANTS[cipaordv1_A41]*exp( CONSTANTS[cipaordv1_B41]*STATES[cipaordv1_V])*STATES[cipaordv1_O]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q41]))/10.0000))) - ( CONSTANTS[cipaordv1_A52]*exp( CONSTANTS[cipaordv1_B52]*STATES[cipaordv1_V])*STATES[cipaordv1_C2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q52]))/10.0000) -  CONSTANTS[cipaordv1_A62]*exp( CONSTANTS[cipaordv1_B62]*STATES[cipaordv1_V])*STATES[cipaordv1_IC2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q62]))/10.0000));
RATES[cipaordv1_O] = (( CONSTANTS[cipaordv1_A31]*exp( CONSTANTS[cipaordv1_B31]*STATES[cipaordv1_V])*STATES[cipaordv1_C2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q31]))/10.0000) -  CONSTANTS[cipaordv1_A41]*exp( CONSTANTS[cipaordv1_B41]*STATES[cipaordv1_V])*STATES[cipaordv1_O]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q41]))/10.0000)) - ( CONSTANTS[cipaordv1_A53]*exp( CONSTANTS[cipaordv1_B53]*STATES[cipaordv1_V])*STATES[cipaordv1_O]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q53]))/10.0000) -  CONSTANTS[cipaordv1_A63]*exp( CONSTANTS[cipaordv1_B63]*STATES[cipaordv1_V])*STATES[cipaordv1_IO]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q63]))/10.0000))) - ( (( CONSTANTS[cipaordv1_Kmax]*CONSTANTS[cipaordv1_Ku]*pow( STATES[cipaordv1_D],CONSTANTS[cipaordv1_n]))/(pow( STATES[cipaordv1_D],CONSTANTS[cipaordv1_n])+CONSTANTS[cipaordv1_halfmax]))*STATES[cipaordv1_O] -  CONSTANTS[cipaordv1_Ku]*STATES[cipaordv1_Obound]);
RATES[cipaordv1_IO] = ((( CONSTANTS[cipaordv1_A3]*exp( CONSTANTS[cipaordv1_B3]*STATES[cipaordv1_V])*STATES[cipaordv1_IC2]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q3]))/10.0000) -  CONSTANTS[cipaordv1_A4]*exp( CONSTANTS[cipaordv1_B4]*STATES[cipaordv1_V])*STATES[cipaordv1_IO]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q4]))/10.0000))+ CONSTANTS[cipaordv1_A53]*exp( CONSTANTS[cipaordv1_B53]*STATES[cipaordv1_V])*STATES[cipaordv1_O]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q53]))/10.0000)) -  CONSTANTS[cipaordv1_A63]*exp( CONSTANTS[cipaordv1_B63]*STATES[cipaordv1_V])*STATES[cipaordv1_IO]*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q63]))/10.0000)) - ( (( CONSTANTS[cipaordv1_Kmax]*CONSTANTS[cipaordv1_Ku]*pow( STATES[cipaordv1_D],CONSTANTS[cipaordv1_n]))/(pow( STATES[cipaordv1_D],CONSTANTS[cipaordv1_n])+CONSTANTS[cipaordv1_halfmax]))*STATES[cipaordv1_IO] -  (( CONSTANTS[cipaordv1_Ku]*CONSTANTS[cipaordv1_A53]*exp( CONSTANTS[cipaordv1_B53]*STATES[cipaordv1_V])*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q53]))/10.0000))/( CONSTANTS[cipaordv1_A63]*exp( CONSTANTS[cipaordv1_B63]*STATES[cipaordv1_V])*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q63]))/10.0000)))*STATES[cipaordv1_IObound]);
RATES[cipaordv1_IObound] = (( (( CONSTANTS[cipaordv1_Kmax]*CONSTANTS[cipaordv1_Ku]*pow( STATES[cipaordv1_D],CONSTANTS[cipaordv1_n]))/(pow( STATES[cipaordv1_D],CONSTANTS[cipaordv1_n])+CONSTANTS[cipaordv1_halfmax]))*STATES[cipaordv1_IO] -  (( CONSTANTS[cipaordv1_Ku]*CONSTANTS[cipaordv1_A53]*exp( CONSTANTS[cipaordv1_B53]*STATES[cipaordv1_V])*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q53]))/10.0000))/( CONSTANTS[cipaordv1_A63]*exp( CONSTANTS[cipaordv1_B63]*STATES[cipaordv1_V])*exp(( (CONSTANTS[cipaordv1_Temp] - 20.0000)*log(CONSTANTS[cipaordv1_q63]))/10.0000)))*STATES[cipaordv1_IObound])+ (CONSTANTS[cipaordv1_Kt]/(1.00000+exp(- (STATES[cipaordv1_V] - CONSTANTS[cipaordv1_Vhalf])/6.78900)))*STATES[cipaordv1_Cbound]) -  CONSTANTS[cipaordv1_Kt]*STATES[cipaordv1_IObound];
RATES[cipaordv1_Obound] = (( (( CONSTANTS[cipaordv1_Kmax]*CONSTANTS[cipaordv1_Ku]*pow( STATES[cipaordv1_D],CONSTANTS[cipaordv1_n]))/(pow( STATES[cipaordv1_D],CONSTANTS[cipaordv1_n])+CONSTANTS[cipaordv1_halfmax]))*STATES[cipaordv1_O] -  CONSTANTS[cipaordv1_Ku]*STATES[cipaordv1_Obound])+ (CONSTANTS[cipaordv1_Kt]/(1.00000+exp(- (STATES[cipaordv1_V] - CONSTANTS[cipaordv1_Vhalf])/6.78900)))*STATES[cipaordv1_Cbound]) -  CONSTANTS[cipaordv1_Kt]*STATES[cipaordv1_Obound];
RATES[cipaordv1_Cbound] = - ( (CONSTANTS[cipaordv1_Kt]/(1.00000+exp(- (STATES[cipaordv1_V] - CONSTANTS[cipaordv1_Vhalf])/6.78900)))*STATES[cipaordv1_Cbound] -  CONSTANTS[cipaordv1_Kt]*STATES[cipaordv1_Obound]) - ( (CONSTANTS[cipaordv1_Kt]/(1.00000+exp(- (STATES[cipaordv1_V] - CONSTANTS[cipaordv1_Vhalf])/6.78900)))*STATES[cipaordv1_Cbound] -  CONSTANTS[cipaordv1_Kt]*STATES[cipaordv1_IObound]);
RATES[cipaordv1_hL] = (ALGEBRAIC[cipaordv1_hLss] - STATES[cipaordv1_hL])/CONSTANTS[cipaordv1_thL];
RATES[cipaordv1_hLp] = (ALGEBRAIC[cipaordv1_hLssp] - STATES[cipaordv1_hLp])/CONSTANTS[cipaordv1_thLp];
RATES[cipaordv1_m] = (ALGEBRAIC[cipaordv1_mss] - STATES[cipaordv1_m])/ALGEBRAIC[cipaordv1_tm];
RATES[cipaordv1_hf] = (ALGEBRAIC[cipaordv1_hss] - STATES[cipaordv1_hf])/ALGEBRAIC[cipaordv1_thf];
RATES[cipaordv1_hs] = (ALGEBRAIC[cipaordv1_hss] - STATES[cipaordv1_hs])/ALGEBRAIC[cipaordv1_ths];
RATES[cipaordv1_a] = (ALGEBRAIC[cipaordv1_ass] - STATES[cipaordv1_a])/ALGEBRAIC[cipaordv1_ta];
RATES[cipaordv1_d] = (ALGEBRAIC[cipaordv1_dss] - STATES[cipaordv1_d])/ALGEBRAIC[cipaordv1_td];
RATES[cipaordv1_ff] = (ALGEBRAIC[cipaordv1_fss] - STATES[cipaordv1_ff])/ALGEBRAIC[cipaordv1_tff];
RATES[cipaordv1_fs] = (ALGEBRAIC[cipaordv1_fss] - STATES[cipaordv1_fs])/ALGEBRAIC[cipaordv1_tfs];
RATES[cipaordv1_jca] = (ALGEBRAIC[cipaordv1_fcass] - STATES[cipaordv1_jca])/CONSTANTS[cipaordv1_tjca];
RATES[cipaordv1_nca] =  ALGEBRAIC[cipaordv1_anca]*CONSTANTS[cipaordv1_k2n] -  STATES[cipaordv1_nca]*ALGEBRAIC[cipaordv1_km2n];
RATES[cipaordv1_xs1] = (ALGEBRAIC[cipaordv1_xs1ss] - STATES[cipaordv1_xs1])/ALGEBRAIC[cipaordv1_txs1];
RATES[cipaordv1_xk1] = (ALGEBRAIC[cipaordv1_xk1ss] - STATES[cipaordv1_xk1])/ALGEBRAIC[cipaordv1_txk1];
RATES[cipaordv1_CaMKt] =  CONSTANTS[cipaordv1_aCaMK]*ALGEBRAIC[cipaordv1_CaMKb]*(ALGEBRAIC[cipaordv1_CaMKb]+STATES[cipaordv1_CaMKt]) -  CONSTANTS[cipaordv1_bCaMK]*STATES[cipaordv1_CaMKt];
RATES[cipaordv1_j] = (ALGEBRAIC[cipaordv1_jss] - STATES[cipaordv1_j])/ALGEBRAIC[cipaordv1_tj];
RATES[cipaordv1_ap] = (ALGEBRAIC[cipaordv1_assp] - STATES[cipaordv1_ap])/ALGEBRAIC[cipaordv1_ta];
RATES[cipaordv1_fcaf] = (ALGEBRAIC[cipaordv1_fcass] - STATES[cipaordv1_fcaf])/ALGEBRAIC[cipaordv1_tfcaf];
RATES[cipaordv1_fcas] = (ALGEBRAIC[cipaordv1_fcass] - STATES[cipaordv1_fcas])/ALGEBRAIC[cipaordv1_tfcas];
RATES[cipaordv1_ffp] = (ALGEBRAIC[cipaordv1_fss] - STATES[cipaordv1_ffp])/ALGEBRAIC[cipaordv1_tffp];
RATES[cipaordv1_xs2] = (ALGEBRAIC[cipaordv1_xs2ss] - STATES[cipaordv1_xs2])/ALGEBRAIC[cipaordv1_txs2];
RATES[cipaordv1_hsp] = (ALGEBRAIC[cipaordv1_hssp] - STATES[cipaordv1_hsp])/ALGEBRAIC[cipaordv1_thsp];
RATES[cipaordv1_jp] = (ALGEBRAIC[cipaordv1_jss] - STATES[cipaordv1_jp])/ALGEBRAIC[cipaordv1_tjp];
RATES[cipaordv1_mL] = (ALGEBRAIC[cipaordv1_mLss] - STATES[cipaordv1_mL])/ALGEBRAIC[cipaordv1_tmL];
RATES[cipaordv1_fcafp] = (ALGEBRAIC[cipaordv1_fcass] - STATES[cipaordv1_fcafp])/ALGEBRAIC[cipaordv1_tfcafp];
RATES[cipaordv1_iF] = (ALGEBRAIC[cipaordv1_iss] - STATES[cipaordv1_iF])/ALGEBRAIC[cipaordv1_tiF];
RATES[cipaordv1_iS] = (ALGEBRAIC[cipaordv1_iss] - STATES[cipaordv1_iS])/ALGEBRAIC[cipaordv1_tiS];
RATES[cipaordv1_iFp] = (ALGEBRAIC[cipaordv1_iss] - STATES[cipaordv1_iFp])/ALGEBRAIC[cipaordv1_tiFp];
RATES[cipaordv1_iSp] = (ALGEBRAIC[cipaordv1_iss] - STATES[cipaordv1_iSp])/ALGEBRAIC[cipaordv1_tiSp];
RATES[cipaordv1_Jrelnp] = (ALGEBRAIC[cipaordv1_Jrel_inf] - STATES[cipaordv1_Jrelnp])/ALGEBRAIC[cipaordv1_tau_rel];
RATES[cipaordv1_Jrelp] = (ALGEBRAIC[cipaordv1_Jrel_infp] - STATES[cipaordv1_Jrelp])/ALGEBRAIC[cipaordv1_tau_relp];
RATES[cipaordv1_ki] = ( - ((ALGEBRAIC[cipaordv1_Ito]+ALGEBRAIC[cipaordv1_IKr]+ALGEBRAIC[cipaordv1_IKs]+ALGEBRAIC[cipaordv1_IK1]+ALGEBRAIC[cipaordv1_IKb]+ALGEBRAIC[cipaordv1_Istim]) -  2.00000*ALGEBRAIC[cipaordv1_INaK])*CONSTANTS[cipaordv1_cm]*CONSTANTS[cipaordv1_Acap])/( CONSTANTS[cipaordv1_F]*CONSTANTS[cipaordv1_vmyo])+( ALGEBRAIC[cipaordv1_JdiffK]*CONSTANTS[cipaordv1_vss])/CONSTANTS[cipaordv1_vmyo];
RATES[cipaordv1_kss] = ( - ALGEBRAIC[cipaordv1_ICaK]*CONSTANTS[cipaordv1_cm]*CONSTANTS[cipaordv1_Acap])/( CONSTANTS[cipaordv1_F]*CONSTANTS[cipaordv1_vss]) - ALGEBRAIC[cipaordv1_JdiffK];
RATES[cipaordv1_nai] = ( - (ALGEBRAIC[cipaordv1_INa]+ALGEBRAIC[cipaordv1_INaL]+ 3.00000*ALGEBRAIC[cipaordv1_INaCa_i]+ 3.00000*ALGEBRAIC[cipaordv1_INaK]+ALGEBRAIC[cipaordv1_INab])*CONSTANTS[cipaordv1_Acap]*CONSTANTS[cipaordv1_cm])/( CONSTANTS[cipaordv1_F]*CONSTANTS[cipaordv1_vmyo])+( ALGEBRAIC[cipaordv1_JdiffNa]*CONSTANTS[cipaordv1_vss])/CONSTANTS[cipaordv1_vmyo];
RATES[cipaordv1_nass] = ( - (ALGEBRAIC[cipaordv1_ICaNa]+ 3.00000*ALGEBRAIC[cipaordv1_INaCa_ss])*CONSTANTS[cipaordv1_cm]*CONSTANTS[cipaordv1_Acap])/( CONSTANTS[cipaordv1_F]*CONSTANTS[cipaordv1_vss]) - ALGEBRAIC[cipaordv1_JdiffNa];
RATES[cipaordv1_V] = - (ALGEBRAIC[cipaordv1_INa]+ALGEBRAIC[cipaordv1_INaL]+ALGEBRAIC[cipaordv1_Ito]+ALGEBRAIC[cipaordv1_ICaL]+ALGEBRAIC[cipaordv1_ICaNa]+ALGEBRAIC[cipaordv1_ICaK]+ALGEBRAIC[cipaordv1_IKr]+ALGEBRAIC[cipaordv1_IKs]+ALGEBRAIC[cipaordv1_IK1]+ALGEBRAIC[cipaordv1_INaCa_i]+ALGEBRAIC[cipaordv1_INaCa_ss]+ALGEBRAIC[cipaordv1_INaK]+ALGEBRAIC[cipaordv1_INab]+ALGEBRAIC[cipaordv1_IKb]+ALGEBRAIC[cipaordv1_IpCa]+ALGEBRAIC[cipaordv1_ICab]+ALGEBRAIC[cipaordv1_Istim]);
RATES[cipaordv1_cass] =  ALGEBRAIC[cipaordv1_Bcass]*((( - (ALGEBRAIC[cipaordv1_ICaL] -  2.00000*ALGEBRAIC[cipaordv1_INaCa_ss])*CONSTANTS[cipaordv1_cm]*CONSTANTS[cipaordv1_Acap])/( 2.00000*CONSTANTS[cipaordv1_F]*CONSTANTS[cipaordv1_vss])+( ALGEBRAIC[cipaordv1_Jrel]*CONSTANTS[cipaordv1_vjsr])/CONSTANTS[cipaordv1_vss]) - ALGEBRAIC[cipaordv1_Jdiff]);
// new for coupling
RATES[cipaordv1_ca_trpn] = CONSTANTS[cipaordv1_trpnmax] * RATES[land2016_TRPN];
RATES[cipaordv1_cai] =  ALGEBRAIC[cipaordv1_Bcai]*((( - ((ALGEBRAIC[cipaordv1_IpCa]+ALGEBRAIC[cipaordv1_ICab]) -  2.00000*ALGEBRAIC[cipaordv1_INaCa_i])*CONSTANTS[cipaordv1_cm]*CONSTANTS[cipaordv1_Acap])/( 2.00000*CONSTANTS[cipaordv1_F]*CONSTANTS[cipaordv1_vmyo]) - ( ALGEBRAIC[cipaordv1_Jup]*CONSTANTS[cipaordv1_vnsr])/CONSTANTS[cipaordv1_vmyo])+( ALGEBRAIC[cipaordv1_Jdiff]*CONSTANTS[cipaordv1_vss])/CONSTANTS[cipaordv1_vmyo] - RATES[cipaordv1_ca_trpn]); //modified
RATES[cipaordv1_cansr] = ALGEBRAIC[cipaordv1_Jup] - ( ALGEBRAIC[cipaordv1_Jtr]*CONSTANTS[cipaordv1_vjsr])/CONSTANTS[cipaordv1_vnsr];
RATES[cipaordv1_cajsr] =  ALGEBRAIC[cipaordv1_Bcajsr]*(ALGEBRAIC[cipaordv1_Jtr] - ALGEBRAIC[cipaordv1_Jrel]);

}

void cipaordv1_land2016::solveEuler(double dt){
  // std::cout <<"Euler cipaordv1_land2016...\n";
  for(int i=0;i<states_size;i++){
    STATES[i] = STATES[i] + RATES[i] * dt;
  }
}

void cipaordv1_land2016::solveRK4(double TIME, double dt){
	unsigned short idx;
	double k1[49],k2[49],k3[49],k4[49];
	double states_temp[49];

	computeRates(TIME, CONSTANTS, RATES, STATES, ALGEBRAIC );
	for(idx = 0; idx < states_size; idx++){
		k1[idx] = dt * RATES[idx];
		states_temp[idx] = STATES[idx] + k1[idx]*0.5;
	}
	computeRates(TIME+(dt*0.5), CONSTANTS, RATES, states_temp, ALGEBRAIC );
	for(idx = 0; idx < states_size; idx++){
		k2[idx] = dt * RATES[idx];
		states_temp[idx] = STATES[idx] + k2[idx]*0.5;
	}
	computeRates(TIME+(dt*0.5), CONSTANTS, RATES, states_temp, ALGEBRAIC );
	for(idx = 0; idx < states_size; idx++){
		k3[idx] = dt * RATES[idx];
		states_temp[idx] = STATES[idx] + k3[idx];
	}
	computeRates(TIME+dt, CONSTANTS, RATES, states_temp, ALGEBRAIC );
	for(idx = 0; idx < states_size; idx++){
		k4[idx] = dt * RATES[idx];
		STATES[idx] += (k1[idx]/6) + (k2[idx]/3) + (k3[idx]/3) + (k4[idx]/6) ;
	}
}

double cipaordv1_land2016::set_time_step(double TIME,
                                              double time_point,
                                              double min_time_step,
                                              double max_time_step,
                                              double min_dV,
                                              double max_dV,
                                              double* CONSTANTS,
                                              double* RATES,
                                              double* STATES,
                                              double* ALGEBRAIC) {
// double time_step = 0.005;
 double time_step = min_time_step;
 if (TIME <= time_point || (TIME - floor(TIME / CONSTANTS[cipaordv1_BCL]) * CONSTANTS[cipaordv1_BCL]) <= time_point) {
    //printf("TIME <= time_point ms\n");
    return time_step;
    //printf("TIME = %E, dV = %E, time_step = %E\n",TIME, RATES[cipaordv1_V] * time_step, time_step);
  }
  else {
    //printf("TIME > time_point ms\n");
    if (std::abs(RATES[cipaordv1_V] * time_step) <= min_dV) {//Slow changes in V
        //printf("dV/dt <= 0.2\n");
        time_step = std::abs(max_dV / RATES[cipaordv1_V]);
        //Make sure time_step is between min time step and max_time_step
        if (time_step < min_time_step) {
            time_step = min_time_step;
        }
        else if (time_step > max_time_step) {
            time_step = max_time_step;
        }
        //printf("TIME = %E, dV = %E, time_step = %E\n",TIME, RATES[cipaordv1_V] * time_step, time_step);
    }
    else if (std::abs(RATES[cipaordv1_V] * time_step) >= max_dV) {//Fast changes in V
        //printf("dV/dt >= 0.8\n");
        time_step = std::abs(min_dV / RATES[cipaordv1_V]);
        //Make sure time_step is not less than 0.005
        if (time_step < min_time_step) {
            time_step = min_time_step;
        }
        //printf("TIME = %E, dV = %E, time_step = %E\n",TIME, RATES[cipaordv1_V] * time_step, time_step);
    } else {
        time_step = min_time_step;
    }
    return time_step;
  }
}

void cipaordv1_land2016::___applyDrugEffect(double conc, const double *hill)
{
CONSTANTS[cipaordv1_GK1] = CONSTANTS[cipaordv1_GK1] * ((hill[2] > 10E-14 && hill[3] > 10E-14) ? 1./(1.+pow(conc/hill[2],hill[3])) : 1.);
CONSTANTS[cipaordv1_GKs] = CONSTANTS[cipaordv1_GKs] * ((hill[4] > 10E-14 && hill[5] > 10E-14) ? 1./(1.+pow(conc/hill[4],hill[5])) : 1.);
CONSTANTS[cipaordv1_GNaL] = CONSTANTS[cipaordv1_GNaL] * ((hill[8] > 10E-14 && hill[9] > 10E-14) ? 1./(1.+pow(conc/hill[8],hill[9])) : 1.);
CONSTANTS[cipaordv1_GNa] = CONSTANTS[cipaordv1_GNa] * ((hill[6] > 10E-14 && hill[7] > 10E-14) ? 1./(1.+pow(conc/hill[6],hill[7])) : 1.);
CONSTANTS[cipaordv1_Gto] = CONSTANTS[cipaordv1_Gto] * ((hill[10] > 10E-14 && hill[11] > 10E-14) ? 1./(1.+pow(conc/hill[10],hill[11])) : 1.);
CONSTANTS[cipaordv1_PCa] = CONSTANTS[cipaordv1_PCa] * ( (hill[0] > 10E-14 && hill[1] > 10E-14) ? 1./(1.+pow(conc/hill[0],hill[1])) : 1.);
}

void cipaordv1_land2016::___applyHERGBinding(double conc, const double *herg)
{
if(conc > 10E-14){
CONSTANTS[cipaordv1_Kmax] = herg[0];
CONSTANTS[cipaordv1_Ku] = herg[1];
CONSTANTS[cipaordv1_n] = herg[2];
CONSTANTS[cipaordv1_halfmax] = herg[3];
CONSTANTS[cipaordv1_Vhalf] = herg[4];
CONSTANTS[cipaordv1_cnc] = conc;
CONSTANTS[cipaordv1_Kt] = 0.000035;
STATES[cipaordv1_D] = conc;
}
}

