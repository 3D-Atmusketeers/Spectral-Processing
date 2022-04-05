/*Rotation*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "input.h"
#include "opac.h"
#include "atmos.h"
#include "constant.h"
#include "include.h"
#include "nrutil.h"
// C HARADA -- update for 2stream //
#include "two_stream.h"


/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;
extern struct Opac opac;

/* --- Function prototypes --------------------------------------- */

void Locate(int n, double *array, double value, int *ilow);
double Planck(double T, double lambda);
double lint2D(double x1, double x2, double y1, double y2, double z1,
              double z2, double z3, double z4, double x, double y);
double lint3D(double x1, double x2, double y1, double y2, double z1,
              double z2, double f1, double f2, double f3, double f4,
              double f5, double f6, double f7, double f8, double x,
              double y, double z);
void Angles3d(double ds[], double theta[], double dtheta[], double lat);
double Radius(double R_pl, double ds[]);
double lint(double xa, double ya, double xb, double yb, double x);

// C HARADA -- update for 2stream //
double two_stream(int NLAYER, int kmin, double *w0_array, double *g0_array, \
                  const double *temperature_array, const double *tau_array, \
                  double NU, double NU_BIN, double incident_frac, double *dtau_array);

/* ------- begin ---------------- RT_Emit -------------------- */

/* Rest frame or Doppler shifted emission spectra for the disk (not including the limb) with clouds turned on or off */

int RT_Emit_3D(double PHASE)
{
    double ***tau_tr_east, ***tau_tr_west, **theta, **dtheta, ***tau_em, ***dtau_em, ***phi_lon_solid, ***theta_lat_solid, ***temperature_3d,
    ***aero_kappa_pre_qext_1, ***aero_tau_pre_qext_1,
    ***aero_kappa_pre_qext_2, ***aero_tau_pre_qext_2,
    ***aero_kappa_pre_qext_3, ***aero_tau_pre_qext_3,
    ***aero_kappa_pre_qext_4, ***aero_tau_pre_qext_4,
    ***aero_kappa_pre_qext_5, ***aero_tau_pre_qext_5,
    ***aero_kappa_pre_qext_6, ***aero_tau_pre_qext_6,
    ***aero_kappa_pre_qext_7, ***aero_tau_pre_qext_7,
    ***aero_kappa_pre_qext_8, ***aero_tau_pre_qext_8,
    ***aero_kappa_pre_qext_9, ***aero_tau_pre_qext_9,
    ***aero_kappa_pre_qext_10, ***aero_tau_pre_qext_10,
    ***aero_kappa_pre_qext_11, ***aero_tau_pre_qext_11,
    ***aero_kappa_pre_qext_12, ***aero_tau_pre_qext_12,
    ***aero_kappa_pre_qext_13, ***aero_tau_pre_qext_13,
    ***kappa_nu_array, ***pressure_array;
    double **intensity, **bad_interp_array, *flux_st, *flux_pl, *flux_tr, *ds, ***dl, **phi,
    *phi_plus_e, *phi_plus_w, *phi_minus_e, *phi_minus_w, **dphi, *theta_lon, *theta_lat, *phi_lon;
    double R, test, a, b, *lat_rad, kappa_nu_plus_e, kappa_nu_plus_w,
    kappa_nu_minus_e, kappa_nu_minus_w, t_lon_plus_e, t_lon_plus_w,
    t_lon_minus_e, t_lon_minus_w, p_lon_plus_e, p_lon_plus_w,
    p_lon_minus_e, p_lon_minus_w, temperature, pressure, kappa_nu, *lon_rad,
    aero_kappa_pre_qext_interp_1,
    aero_kappa_pre_qext_interp_2,
    aero_kappa_pre_qext_interp_3,
    aero_kappa_pre_qext_interp_4,
    aero_kappa_pre_qext_interp_5,
    aero_kappa_pre_qext_interp_6,
    aero_kappa_pre_qext_interp_7,
    aero_kappa_pre_qext_interp_8,
    aero_kappa_pre_qext_interp_9,
    aero_kappa_pre_qext_interp_10,
    aero_kappa_pre_qext_interp_11,
    aero_kappa_pre_qext_interp_12,
    aero_kappa_pre_qext_interp_13,
    aero_kappa_1,
    aero_kappa_2,
    aero_kappa_3,
    aero_kappa_4,
    aero_kappa_5,
    aero_kappa_6,
    aero_kappa_7,
    aero_kappa_8,
    aero_kappa_9,
    aero_kappa_10,
    aero_kappa_11,
    aero_kappa_12,
    aero_kappa_13;
    double total_cloud_kappa;
    double cloud_param;
    double **I_top, *I_bot, **dkappa_nu;
    int i, j, k, l, m, n, o, c, g, h, ii;
    double dphid, thetad, dthetad;
    FILE *file;
    double solid;
    double average, running_sum, num_points, good_intensity;
    double u_vel, v_vel, w_vel, v_los, delta_lam, omega;
    double weight_1, weight_2, weight_3, weight_4, weight_5, weight_6, weight_7;
    double weight_8, weight_9, weight_10, weight_11, weight_12, weight_13;
    double temp_value;

    double pressure_array_for_scattering_data_in_pascals[50] = {
    10.02, 3319.25, 5768.19, 7603.96, 10023.97, 10023.97, 13214.17, 17419.66, 22963.59, 22963.59, \
    30271.9, 39906.14, 39906.14, 52606.54, 52606.54, 69348.93, 91419.69, 91419.69, 120514.62, \
    120514.62, 158869.22, 158869.22, 209430.41, 209430.41, 276083.05, 276083.05, 363948.35, \
    479777.36, 479777.36, 632469.74, 632469.74, 833757.5, 833757.5, 1099106.44, 1099106.44, \
    1448904.47, 1448904.47, 1910028.08, 1910028.08, 2517907.38, 2517907.38, 3319248.37, \
    3319248.37, 4375621.54, 4375621.54, 5768192.58, 5768192.58, 7603958.73, 7803958.73, \
    8003958.73};

    double wavelengths_in_microns[50] = {0.3, 0.35, 0.42, 0.49, 0.57, 0.68, 0.79, 0.93, 1.1, 1.29, 1.52, 1.79, 2.1, \
                                         2.47, 2.91, 3.42, 4.02, 4.73, 5.57, 6.55, 7.7, 9.05, 10.65, 12.53, 14.73, \
                                         17.33, 20.38, 23.97, 28.19, 33.16, 39.0, 45.87, 53.95, 63.45, 74.63, \
                                         87.78, 103.24, 121.43, 142.82, 167.97, 197.56, 232.37, 273.3, 321.44, \
                                         378.06, 444.66, 522.99, 615.12, 723.47, 850.92};

    int x=0, y=0, num_wavelength_points=0, num_pressure_points=0;
    double input_val=0;
    num_wavelength_points = 50;
    num_pressure_points = 50;

    const char KCl_wav_gg_file[]    = "SCATTERING_DATA/KCl_wav_gg.txt";
    const char KCl_wav_pi0_file[] = "SCATTERING_DATA/KCl_wav_pi0.txt";
    const char KCl_wav_qext_file[] = "SCATTERING_DATA/KCl_wav_qext.txt";

    const char ZnS_wav_gg_file[]    = "SCATTERING_DATA/ZnS_wav_gg.txt";
    const char ZnS_wav_pi0_file[] = "SCATTERING_DATA/ZnS_wav_pi0.txt";
    const char ZnS_wav_qext_file[] = "SCATTERING_DATA/ZnS_wav_qext.txt";

    const char Na2S_wav_gg_file[]    = "SCATTERING_DATA/Na2S_wav_gg.txt";
    const char Na2S_wav_pi0_file[] = "SCATTERING_DATA/Na2S_wav_pi0.txt";
    const char Na2S_wav_qext_file[] = "SCATTERING_DATA/Na2S_wav_qext.txt";

    const char MnS_wav_gg_file[]    = "SCATTERING_DATA/MnS_wav_gg.txt";
    const char MnS_wav_pi0_file[] = "SCATTERING_DATA/MnS_wav_pi0.txt";
    const char MnS_wav_qext_file[] = "SCATTERING_DATA/MnS_wav_qext.txt";

    const char Cr_wav_gg_file[]    = "SCATTERING_DATA/Cr_wav_gg.txt";
    const char Cr_wav_pi0_file[] = "SCATTERING_DATA/Cr_wav_pi0.txt";
    const char Cr_wav_qext_file[] = "SCATTERING_DATA/Cr_wav_qext.txt";

    const char SiO2_wav_gg_file[]    = "SCATTERING_DATA/SiO2_wav_gg.txt";
    const char SiO2_wav_pi0_file[] = "SCATTERING_DATA/SiO2_wav_pi0.txt";
    const char SiO2_wav_qext_file[] = "SCATTERING_DATA/SiO2_wav_qext.txt";

    const char Mg2SiO4_wav_gg_file[]    = "SCATTERING_DATA/Mg2SiO4_wav_gg.txt";
    const char Mg2SiO4_wav_pi0_file[] = "SCATTERING_DATA/Mg2SiO4_wav_pi0.txt";
    const char Mg2SiO4_wav_qext_file[] = "SCATTERING_DATA/Mg2SiO4_wav_qext.txt";

    const char VO_wav_gg_file[]    = "SCATTERING_DATA/VO_wav_gg.txt";
    const char VO_wav_pi0_file[] = "SCATTERING_DATA/VO_wav_pi0.txt";
    const char VO_wav_qext_file[] = "SCATTERING_DATA/VO_wav_qext.txt";

    const char Ni_wav_gg_file[]    = "SCATTERING_DATA/Ni_wav_gg.txt";
    const char Ni_wav_pi0_file[] = "SCATTERING_DATA/Ni_wav_pi0.txt";
    const char Ni_wav_qext_file[] = "SCATTERING_DATA/Ni_wav_qext.txt";

    const char Fe_wav_gg_file[]    = "SCATTERING_DATA/Fe_wav_gg.txt";
    const char Fe_wav_pi0_file[] = "SCATTERING_DATA/Fe_wav_pi0.txt";
    const char Fe_wav_qext_file[] = "SCATTERING_DATA/Fe_wav_qext.txt";

    const char CaSiO4_wav_gg_file[]    = "SCATTERING_DATA/CaSiO4_wav_gg.txt";
    const char CaSiO4_wav_pi0_file[] = "SCATTERING_DATA/CaSiO4_wav_pi0.txt";
    const char CaSiO4_wav_qext_file[] = "SCATTERING_DATA/CaSiO4_wav_qext.txt";

    const char CaTiO3_wav_gg_file[]    = "SCATTERING_DATA/CaTiO3_wav_gg.txt";
    const char CaTiO3_wav_pi0_file[] = "SCATTERING_DATA/CaTiO3_wav_pi0.txt";
    const char CaTiO3_wav_qext_file[] = "SCATTERING_DATA/CaTiO3_wav_qext.txt";

    const char Al2O3_wav_gg_file[]    = "SCATTERING_DATA/Al2O3_wav_gg.txt";
    const char Al2O3_wav_pi0_file[] = "SCATTERING_DATA/Al2O3_wav_pi0.txt";
    const char Al2O3_wav_qext_file[] = "SCATTERING_DATA/Al2O3_wav_qext.txt";

    FILE *input_KCl_wav_gg_file;
    FILE *input_KCl_wav_pi0_file;
    FILE *input_KCl_wav_qext_file;

    FILE *input_ZnS_wav_gg_file;
    FILE *input_ZnS_wav_pi0_file;
    FILE *input_ZnS_wav_qext_file;

    FILE *input_Na2S_wav_gg_file;
    FILE *input_Na2S_wav_pi0_file;
    FILE *input_Na2S_wav_qext_file;

    FILE *input_MnS_wav_gg_file;
    FILE *input_MnS_wav_pi0_file;
    FILE *input_MnS_wav_qext_file;

    FILE *input_Cr_wav_gg_file;
    FILE *input_Cr_wav_pi0_file;
    FILE *input_Cr_wav_qext_file;

    FILE *input_SiO2_wav_gg_file;
    FILE *input_SiO2_wav_pi0_file;
    FILE *input_SiO2_wav_qext_file;

    FILE *input_Mg2SiO4_wav_gg_file;
    FILE *input_Mg2SiO4_wav_pi0_file;
    FILE *input_Mg2SiO4_wav_qext_file;

    FILE *input_VO_wav_gg_file;
    FILE *input_VO_wav_pi0_file;
    FILE *input_VO_wav_qext_file;

    FILE *input_Ni_wav_gg_file;
    FILE *input_Ni_wav_pi0_file;
    FILE *input_Ni_wav_qext_file;

    FILE *input_Fe_wav_gg_file;
    FILE *input_Fe_wav_pi0_file;
    FILE *input_Fe_wav_qext_file;

    FILE *input_CaSiO4_wav_gg_file;
    FILE *input_CaSiO4_wav_pi0_file;
    FILE *input_CaSiO4_wav_qext_file;

    FILE *input_CaTiO3_wav_gg_file;
    FILE *input_CaTiO3_wav_pi0_file;
    FILE *input_CaTiO3_wav_qext_file;

    FILE *input_Al2O3_wav_gg_file;
    FILE *input_Al2O3_wav_pi0_file;
    FILE *input_Al2O3_wav_qext_file;

    input_KCl_wav_gg_file = fopen(KCl_wav_gg_file, "r");
    input_KCl_wav_pi0_file = fopen(KCl_wav_pi0_file, "r");
    input_KCl_wav_qext_file = fopen(KCl_wav_qext_file, "r");

    input_ZnS_wav_gg_file = fopen(ZnS_wav_gg_file, "r");
    input_ZnS_wav_pi0_file = fopen(ZnS_wav_pi0_file, "r");
    input_ZnS_wav_qext_file = fopen(ZnS_wav_qext_file, "r");

    input_Na2S_wav_gg_file = fopen(Na2S_wav_gg_file, "r");
    input_Na2S_wav_pi0_file = fopen(Na2S_wav_pi0_file, "r");
    input_Na2S_wav_qext_file = fopen(Na2S_wav_qext_file, "r");

    input_MnS_wav_gg_file = fopen(MnS_wav_gg_file, "r");
    input_MnS_wav_pi0_file = fopen(MnS_wav_pi0_file, "r");
    input_MnS_wav_qext_file = fopen(MnS_wav_qext_file, "r");

    input_Cr_wav_gg_file = fopen(Cr_wav_gg_file, "r");
    input_Cr_wav_pi0_file = fopen(Cr_wav_pi0_file, "r");
    input_Cr_wav_qext_file = fopen(Cr_wav_qext_file, "r");

    input_SiO2_wav_gg_file = fopen(SiO2_wav_gg_file, "r");
    input_SiO2_wav_pi0_file = fopen(SiO2_wav_pi0_file, "r");
    input_SiO2_wav_qext_file = fopen(SiO2_wav_qext_file, "r");

    input_Mg2SiO4_wav_gg_file = fopen(Mg2SiO4_wav_gg_file, "r");
    input_Mg2SiO4_wav_pi0_file = fopen(Mg2SiO4_wav_pi0_file, "r");
    input_Mg2SiO4_wav_qext_file = fopen(Mg2SiO4_wav_qext_file, "r");

    input_VO_wav_gg_file = fopen(VO_wav_gg_file, "r");
    input_VO_wav_pi0_file = fopen(VO_wav_pi0_file, "r");
    input_VO_wav_qext_file = fopen(VO_wav_qext_file, "r");

    input_Ni_wav_gg_file = fopen(Ni_wav_gg_file, "r");
    input_Ni_wav_pi0_file = fopen(Ni_wav_pi0_file, "r");
    input_Ni_wav_qext_file = fopen(Ni_wav_qext_file, "r");

    input_Fe_wav_gg_file = fopen(Fe_wav_gg_file, "r");
    input_Fe_wav_pi0_file = fopen(Fe_wav_pi0_file, "r");
    input_Fe_wav_qext_file = fopen(Fe_wav_qext_file, "r");

    input_CaSiO4_wav_gg_file = fopen(CaSiO4_wav_gg_file, "r");
    input_CaSiO4_wav_pi0_file = fopen(CaSiO4_wav_pi0_file, "r");
    input_CaSiO4_wav_qext_file = fopen(CaSiO4_wav_qext_file, "r");

    input_CaTiO3_wav_gg_file = fopen(CaTiO3_wav_gg_file, "r");
    input_CaTiO3_wav_pi0_file = fopen(CaTiO3_wav_pi0_file, "r");
    input_CaTiO3_wav_qext_file = fopen(CaTiO3_wav_qext_file, "r");

    input_Al2O3_wav_gg_file = fopen(Al2O3_wav_gg_file, "r");
    input_Al2O3_wav_pi0_file = fopen(Al2O3_wav_pi0_file, "r");
    input_Al2O3_wav_qext_file = fopen(Al2O3_wav_qext_file, "r");

    double KCl_wav_gg[num_pressure_points][num_wavelength_points];
    double KCl_wav_pi0[num_pressure_points][num_wavelength_points];
    double KCl_wav_qext[num_pressure_points][num_wavelength_points];

    double ZnS_wav_gg[num_pressure_points][num_wavelength_points];
    double ZnS_wav_pi0[num_pressure_points][num_wavelength_points];
    double ZnS_wav_qext[num_pressure_points][num_wavelength_points];

    double Na2S_wav_gg[num_pressure_points][num_wavelength_points];
    double Na2S_wav_pi0[num_pressure_points][num_wavelength_points];
    double Na2S_wav_qext[num_pressure_points][num_wavelength_points];

    double MnS_wav_gg[num_pressure_points][num_wavelength_points];
    double MnS_wav_pi0[num_pressure_points][num_wavelength_points];
    double MnS_wav_qext[num_pressure_points][num_wavelength_points];

    double Cr_wav_gg[num_pressure_points][num_wavelength_points];
    double Cr_wav_pi0[num_pressure_points][num_wavelength_points];
    double Cr_wav_qext[num_pressure_points][num_wavelength_points];

    double SiO2_wav_gg[num_pressure_points][num_wavelength_points];
    double SiO2_wav_pi0[num_pressure_points][num_wavelength_points];
    double SiO2_wav_qext[num_pressure_points][num_wavelength_points];

    double Mg2SiO4_wav_gg[num_pressure_points][num_wavelength_points];
    double Mg2SiO4_wav_pi0[num_pressure_points][num_wavelength_points];
    double Mg2SiO4_wav_qext[num_pressure_points][num_wavelength_points];

    double VO_wav_gg[num_pressure_points][num_wavelength_points];
    double VO_wav_pi0[num_pressure_points][num_wavelength_points];
    double VO_wav_qext[num_pressure_points][num_wavelength_points];

    double Ni_wav_gg[num_pressure_points][num_wavelength_points];
    double Ni_wav_pi0[num_pressure_points][num_wavelength_points];
    double Ni_wav_qext[num_pressure_points][num_wavelength_points];

    double Fe_wav_gg[num_pressure_points][num_wavelength_points];
    double Fe_wav_pi0[num_pressure_points][num_wavelength_points];
    double Fe_wav_qext[num_pressure_points][num_wavelength_points];

    double CaSiO4_wav_gg[num_pressure_points][num_wavelength_points];
    double CaSiO4_wav_pi0[num_pressure_points][num_wavelength_points];
    double CaSiO4_wav_qext[num_pressure_points][num_wavelength_points];

    double CaTiO3_wav_gg[num_pressure_points][num_wavelength_points];
    double CaTiO3_wav_pi0[num_pressure_points][num_wavelength_points];
    double CaTiO3_wav_qext[num_pressure_points][num_wavelength_points];

    double Al2O3_wav_gg[num_pressure_points][num_wavelength_points];
    double Al2O3_wav_pi0[num_pressure_points][num_wavelength_points];
    double Al2O3_wav_qext[num_pressure_points][num_wavelength_points];


    for(x = 0; x < num_pressure_points; x++)
    {
        for(y = 0; y < num_wavelength_points; y++)
        {
            fscanf(input_KCl_wav_gg_file, "%le", &input_val);
            KCl_wav_gg[x][y]=input_val;
            fscanf(input_KCl_wav_pi0_file, "%le", &input_val);
            KCl_wav_pi0[x][y]=input_val;
            fscanf(input_KCl_wav_qext_file, "%le", &input_val);
            KCl_wav_qext[x][y]=input_val;


            fscanf(input_ZnS_wav_gg_file, "%le", &input_val);
            ZnS_wav_gg[x][y]=input_val;
            fscanf(input_ZnS_wav_pi0_file, "%le", &input_val);
            ZnS_wav_pi0[x][y]=input_val;
            fscanf(input_ZnS_wav_qext_file, "%le", &input_val);
            ZnS_wav_qext[x][y]=input_val;


            fscanf(input_Na2S_wav_gg_file, "%le", &input_val);
            Na2S_wav_gg[x][y]=input_val;
            fscanf(input_Na2S_wav_pi0_file, "%le", &input_val);
            Na2S_wav_pi0[x][y]=input_val;
            fscanf(input_Na2S_wav_qext_file, "%le", &input_val);
            Na2S_wav_qext[x][y]=input_val;


            fscanf(input_MnS_wav_gg_file, "%le", &input_val);
            MnS_wav_gg[x][y]=input_val;
            fscanf(input_MnS_wav_pi0_file, "%le", &input_val);
            MnS_wav_pi0[x][y]=input_val;
            fscanf(input_MnS_wav_qext_file, "%le", &input_val);
            MnS_wav_qext[x][y]=input_val;


            fscanf(input_Cr_wav_gg_file, "%le", &input_val);
            Cr_wav_gg[x][y]=input_val;
            fscanf(input_Cr_wav_pi0_file, "%le", &input_val);
            Cr_wav_pi0[x][y]=input_val;
            fscanf(input_Cr_wav_qext_file, "%le", &input_val);
            Cr_wav_qext[x][y]=input_val;


            fscanf(input_SiO2_wav_gg_file, "%le", &input_val);
            SiO2_wav_gg[x][y]=input_val;
            fscanf(input_SiO2_wav_pi0_file, "%le", &input_val);
            SiO2_wav_pi0[x][y]=input_val;
            fscanf(input_SiO2_wav_qext_file, "%le", &input_val);
            SiO2_wav_qext[x][y]=input_val;


            fscanf(input_Mg2SiO4_wav_pi0_file, "%le", &input_val);
            Mg2SiO4_wav_pi0[x][y]=input_val;
            fscanf(input_Mg2SiO4_wav_pi0_file, "%le", &input_val);
            Mg2SiO4_wav_pi0[x][y]=input_val;
            fscanf(input_Mg2SiO4_wav_qext_file, "%le", &input_val);
            Mg2SiO4_wav_qext[x][y]=input_val;


            fscanf(input_VO_wav_gg_file, "%le", &input_val);
            VO_wav_gg[x][y]=input_val;
            fscanf(input_VO_wav_pi0_file, "%le", &input_val);
            VO_wav_pi0[x][y]=input_val;
            fscanf(input_VO_wav_qext_file, "%le", &input_val);
            VO_wav_qext[x][y]=input_val;


            fscanf(input_Ni_wav_gg_file, "%le", &input_val);
            Ni_wav_gg[x][y]=input_val;
            fscanf(input_Ni_wav_pi0_file, "%le", &input_val);
            Ni_wav_pi0[x][y]=input_val;
            fscanf(input_Ni_wav_qext_file, "%le", &input_val);
            Ni_wav_qext[x][y]=input_val;


            fscanf(input_Fe_wav_gg_file, "%le", &input_val);
            Fe_wav_gg[x][y]=input_val;
            fscanf(input_Fe_wav_pi0_file, "%le", &input_val);
            Fe_wav_pi0[x][y]=input_val;
            fscanf(input_Fe_wav_qext_file, "%le", &input_val);
            Fe_wav_qext[x][y]=input_val;


            fscanf(input_CaSiO4_wav_gg_file, "%le", &input_val);
            CaSiO4_wav_gg[x][y]=input_val;
            fscanf(input_CaSiO4_wav_pi0_file, "%le", &input_val);
            CaSiO4_wav_pi0[x][y]=input_val;
            fscanf(input_CaSiO4_wav_qext_file, "%le", &input_val);
            CaSiO4_wav_qext[x][y]=input_val;


            fscanf(input_CaTiO3_wav_gg_file, "%le", &input_val);
            CaTiO3_wav_gg[x][y]=input_val;
            fscanf(input_CaTiO3_wav_pi0_file, "%le", &input_val);
            CaTiO3_wav_pi0[x][y]=input_val;
            fscanf(input_CaTiO3_wav_qext_file, "%le", &input_val);
            CaTiO3_wav_qext[x][y]=input_val;


            fscanf(input_Al2O3_wav_gg_file, "%le", &input_val);
            Al2O3_wav_gg[x][y]=input_val;
            fscanf(input_Al2O3_wav_pi0_file, "%le", &input_val);
            Al2O3_wav_pi0[x][y]=input_val;
            fscanf(input_Al2O3_wav_qext_file, "%le", &input_val);
            Al2O3_wav_qext[x][y]=input_val;
        }
    }


    // C HARADA -- update for 2stream //
    int kmin, good_l, good_m, good_val, pressure_index, wavelength_index;
    double incident_frac;
    double ***pi0_tot, ***asym_tot;


    char OUTPUT_FILE[200];
    sprintf(OUTPUT_FILE, "%s%06.2f.dat", OUTPUT_PREFIX, PHASE);
    file = fopen(OUTPUT_FILE, "w");

    /*Allocate memory*/
    tau_em = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        tau_em[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            tau_em[l][m] = malloc(NTAU*sizeof(double));
        }
    }

    dtau_em = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        dtau_em[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            dtau_em[l][m] = malloc(NTAU*sizeof(double));
        }
    }

    temperature_3d = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        temperature_3d[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            temperature_3d[l][m] = malloc(NTAU*sizeof(double));
        }
    }

    phi_lon_solid = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        phi_lon_solid[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            phi_lon_solid[l][m] = malloc(NTAU*sizeof(double));
        }
    }

    theta_lat_solid = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        theta_lat_solid[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            theta_lat_solid[l][m] = malloc(NTAU*sizeof(double));
        }
    }

    dl = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        dl[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            dl[l][m] = malloc(NTAU*sizeof(double));
        }
    }


    kappa_nu_array = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        kappa_nu_array[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            kappa_nu_array[l][m] = malloc(NTAU*sizeof(double));
        }
    }

    pressure_array = malloc(NLAT*sizeof(double));
    for(l=0; l<NLAT; l++)
    {
        pressure_array[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            pressure_array[l][m] = malloc(NTAU*sizeof(double));
        }
    }


    // C HARADA -- update for 2stream //

    /* allocate memory for scattering parameters */

    pi0_tot = malloc(NLAT*sizeof(double)); // total single scattering albedo
    for(l=0; l<NLAT; l++)
    {
        pi0_tot[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            pi0_tot[l][m] = malloc(NTAU*sizeof(double));
        }
    }

    asym_tot = malloc(NLAT*sizeof(double)); // total asymmetry parameter
    for(l=0; l<NLAT; l++)
    {
        asym_tot[l] = malloc(NLON*sizeof(double));
        for(m=0; m<NLON; m++)
        {
            asym_tot[l][m] = malloc(NTAU*sizeof(double));
        }
    }


    /* allocate memory for aero taus and kappas if clouds on */

    if(CLOUDS==1){

        /* MgSiO3 */
        aero_kappa_pre_qext_1 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_1[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_1[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_1 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_1[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_1[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        /* Fe */
        aero_kappa_pre_qext_2 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_2[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_2[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_2 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_2[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_2[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        aero_kappa_pre_qext_3 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_3[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_3[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_3 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_3[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_3[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        aero_kappa_pre_qext_4 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_4[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_4[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_4 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_4[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_4[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        aero_kappa_pre_qext_5 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_5[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_5[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_5 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_5[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_5[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        aero_kappa_pre_qext_6 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_6[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_6[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_6 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_6[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_6[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        aero_kappa_pre_qext_7 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_7[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_7[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_7 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_7[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_7[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        aero_kappa_pre_qext_8 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_8[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_8[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_8 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_8[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_8[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        aero_kappa_pre_qext_9 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_9[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_9[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_9 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_9[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_9[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        aero_kappa_pre_qext_10 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_10[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_10[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_10 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_10[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_10[l][m] = malloc(NTAU*sizeof(double));
            }
        }


        aero_kappa_pre_qext_11 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_11[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_11[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_11 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_11[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_11[l][m] = malloc(NTAU*sizeof(double));
            }
        }

        aero_kappa_pre_qext_12 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_12[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_12[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_12 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_12[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_12[l][m] = malloc(NTAU*sizeof(double));
            }
        }


        aero_kappa_pre_qext_13 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_kappa_pre_qext_13[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_kappa_pre_qext_13[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        aero_tau_pre_qext_13 = malloc(NLAT*sizeof(double));
        for(l=0; l<NLAT; l++)
        {
            aero_tau_pre_qext_13[l] = malloc(NLON*sizeof(double));
            for(m=0; m<NLON; m++)
            {
                aero_tau_pre_qext_13[l][m] = malloc(NTAU*sizeof(double));
            }
        }
        printf("clouds: ON\n");
    }

    else{
        printf("clouds: OFF\n");
    }

    if(DOPPLER==1){
        printf("doppler: ON\n");
    }
    else{
        printf("doppler: OFF\n");
    }

    theta = dmatrix(0, NLAT-1, 0, NLON-1);
    dtheta = dmatrix(0, NLAT-1, 0, NLON-1);
    dphi = dmatrix(0, NLAT-1, 0, NLON-1);
    phi = dmatrix(0, NLAT-1, 0, NLON-1);

    intensity = dmatrix(0, NLAT-1, 0, NLON-1);
    bad_interp_array = dmatrix(0, NLAT-1, 0, NLON-1);
    I_top = dmatrix(0, NLAT-1, 0, NLON-1);

    flux_pl = dvector(0, NLAMBDA-1);
    ds = dvector(0, NTAU-1);

    lat_rad = dvector(0, NLAT-1);
    lon_rad = dvector(0, NLON-1);

    /*   Calculate the angular rotational speed omega */

    omega = 2.0*PI / (P_ROT*24.0*60.0*60.0);

    /*Calculate ds*/

    for(j=NTAU-1; j>=0; j--)
    {
        ds[j] = atmos.alt[j-1] - atmos.alt[j];
    }
    ds[0] = ds[1];

    /* calculate new aerosol taus and kappas, corrected for wavelength (if clouds on) */

    if(CLOUDS==1){
        for(l=0; l<NLAT; l++){
            for(m=0; m<NLON; m++){
                for(j=0; j<NTAU; j++){
                    /* scattering efficiency correction from 5um to 2.3um
                     (PI0, G0, QE calculated from Mie Scattering code of Mischenko, used in Roman Malsky 2022) */
                    aero_tau_pre_qext_1[l][m][j] = atmos.aero_tau_pre_qext_1[l][m][j];
                    aero_tau_pre_qext_2[l][m][j] = atmos.aero_tau_pre_qext_2[l][m][j];
                    aero_tau_pre_qext_3[l][m][j] = atmos.aero_tau_pre_qext_3[l][m][j];
                    aero_tau_pre_qext_4[l][m][j] = atmos.aero_tau_pre_qext_4[l][m][j];
                    aero_tau_pre_qext_5[l][m][j] = atmos.aero_tau_pre_qext_5[l][m][j];
                    aero_tau_pre_qext_6[l][m][j] = atmos.aero_tau_pre_qext_6[l][m][j];
                    aero_tau_pre_qext_7[l][m][j] = atmos.aero_tau_pre_qext_7[l][m][j];
                    aero_tau_pre_qext_8[l][m][j] = atmos.aero_tau_pre_qext_8[l][m][j];
                    aero_tau_pre_qext_9[l][m][j] = atmos.aero_tau_pre_qext_9[l][m][j];
                    aero_tau_pre_qext_10[l][m][j] = atmos.aero_tau_pre_qext_10[l][m][j];
                    aero_tau_pre_qext_11[l][m][j] = atmos.aero_tau_pre_qext_11[l][m][j];
                    aero_tau_pre_qext_12[l][m][j] = atmos.aero_tau_pre_qext_12[l][m][j];
                    aero_tau_pre_qext_13[l][m][j] = atmos.aero_tau_pre_qext_13[l][m][j];
                }
            }
        }

        for(l=0; l<NLAT; l++){
            for(m=0; m<NLON; m++){
                for(j=0; j<NTAU; j++){
                    aero_kappa_pre_qext_1[l][m][j] = aero_tau_pre_qext_1[l][m][j] / ds[j];
                    aero_kappa_pre_qext_2[l][m][j] = aero_tau_pre_qext_2[l][m][j] / ds[j];
                    aero_kappa_pre_qext_3[l][m][j] = aero_tau_pre_qext_3[l][m][j] / ds[j];
                    aero_kappa_pre_qext_4[l][m][j] = aero_tau_pre_qext_4[l][m][j] / ds[j];
                    aero_kappa_pre_qext_5[l][m][j] = aero_tau_pre_qext_5[l][m][j] / ds[j];
                    aero_kappa_pre_qext_6[l][m][j] = aero_tau_pre_qext_6[l][m][j] / ds[j];
                    aero_kappa_pre_qext_7[l][m][j] = aero_tau_pre_qext_7[l][m][j] / ds[j];
                    aero_kappa_pre_qext_8[l][m][j] = aero_tau_pre_qext_8[l][m][j] / ds[j];
                    aero_kappa_pre_qext_9[l][m][j] = aero_tau_pre_qext_9[l][m][j] / ds[j];
                    aero_kappa_pre_qext_10[l][m][j] = aero_tau_pre_qext_10[l][m][j] / ds[j];
                    aero_kappa_pre_qext_11[l][m][j] = aero_tau_pre_qext_11[l][m][j] / ds[j];
                    aero_kappa_pre_qext_12[l][m][j] = aero_tau_pre_qext_12[l][m][j] / ds[j];
                    aero_kappa_pre_qext_13[l][m][j] = aero_tau_pre_qext_13[l][m][j] / ds[j];
                }
            }
        }
    }
    /*Geometry*/

    /*Calculating dl longitude and latitude along line-of-sight*/
    for(l=0;l<NLAT;l++)
    {
        lat_rad[l] = atmos.lat[l] * PI/180.0; /*theta, latitude*/

        for(m=0;m<NLON;m++)
        {
            if(atmos.lon[m]>=90.0-PHASE && atmos.lon[m]<=270.0-PHASE)
            {

                lon_rad[m] = atmos.lon[m] * PI/180.0; /*phi, longitude*/

                b = R_PLANET;
                for(j=NTAU-1; j>=0; j--)
                {
                    b += ds[j];

                    dl[l][m][j] = pow(SQ(b) - SQ(R_PLANET * sin(lat_rad[l])) - SQ(R_PLANET * cos(lat_rad[l]) * cos(-PI/2. + lon_rad[m] + PHASE*PI/180.0)), 0.5);
                    phi_lon_solid[l][m][j] = 90.0 + acos((R_PLANET * cos(-PI/2. + lon_rad[m] + PHASE*PI/180.0) * cos(lat_rad[l]))/(pow(SQ(b)-SQ(R_PLANET * sin(lat_rad[l])), 0.5))) * 180.0/PI;
                    theta_lat_solid[l][m][j] = asin((R_PLANET * sin(lat_rad[l]))/b) * 180.0/PI;
                }
                for(j=0; j<NTAU; j++)
                {
                    if(j!=NTAU-1)
                    {
                        dl[l][m][j] -= dl[l][m][j+1];
                    }
                    else
                    {
                        dl[l][m][j] = dl[l][m][j] - R_PLANET * cos(lat_rad[l]) * sin(lon_rad[m] + PHASE*PI/180.0 - PI/2.);
                    }
                }
            }
        }
    }

    /*Calculating solid angle along NLAT*NLON */
    solid = 0.0;
    for(l=0; l<NLAT; l++)
    {
        for(m=0; m<NLON; m++)
        {
            if(atmos.lon[m]>=90.0-PHASE && atmos.lon[m]<=270.0-PHASE)
            {
                if(l<NLAT-1 && l>0)
                {
                    theta[l][m] = theta_lat_solid[l][m][0] * PI/180.0;
                    dtheta[l][m] = 0.5*(theta_lat_solid[l+1][m][0] - theta_lat_solid[l-1][m][0])* PI/180.0;
                }
                theta[0][m] = theta_lat_solid[0][m][0]* PI/180.0;
                theta[NLAT-1][m] = theta_lat_solid[NLAT-1][m][0]* PI/180.0;

                dtheta[0][m] = ( 0.5*(theta_lat_solid[1][m][0]+theta_lat_solid[0][m][0]) + 90)* PI/180.0;
                dtheta[NLAT-1][m] = (90 - 0.5*(theta_lat_solid[NLAT-1][m][0]+theta_lat_solid[NLAT-2][m][0]))* PI/180.0;

                if(m>0 && m<NLON-1)
                {
                    if(atmos.lon[m]>90.0-PHASE && atmos.lon[m]<270.0-PHASE)
                    {
                        dphi[l][m] = 0.5*(phi_lon_solid[l][m+1][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 90.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m+1][0] - phi_lon_solid[l][m][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 270.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                }

                else if(m==0)
                {
                    if(atmos.lon[m]>90.0-PHASE && atmos.lon[m]<270.0-PHASE)
                    {
                        dphi[l][m] = 0.5*(phi_lon_solid[l][m+1][0] - phi_lon_solid[l][NLON-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 90.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m+1][0] - phi_lon_solid[l][m][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 270.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m][0] - phi_lon_solid[l][NLON-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                }

                else if(m==NLON-1)
                {
                    if(atmos.lon[m]>90.0-PHASE && atmos.lon[m]<270.0-PHASE)
                    {
                        dphi[l][m] = 0.5*(phi_lon_solid[l][0][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 90.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][0][0] - phi_lon_solid[l][m][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                    else if(atmos.lon[m] == 270.0-PHASE)
                    {
                        dphi[l][m] = (phi_lon_solid[l][m][0] - phi_lon_solid[l][m-1][0])* PI/180.0;
                        phi[l][m] = phi_lon_solid[l][m][0]* PI/180.0;
                    }
                }

                // Every once in a while dphi becomes negative and it should not be
                dtheta[l][m] = fabs(dtheta[l][m]);
                if (dphi[l][m] < 0)
                {
                    dphi[l][m] = 0.01;
                }
                theta[l][m]  = fabs(theta[l][m]);
                solid += SQ(cos(theta[l][m]))*cos(phi[l][m]-PI)*dtheta[l][m]*dphi[l][m];
            }
        }
    }
    printf("solid %f\n", solid);

    for(i=0; i<NLAMBDA; i++)
    {
        // Find the nearest wavelength  index //
        wavelength_index = 0;
        if (atmos.lambda[i] * 1e6 > wavelengths_in_microns[num_wavelength_points - 1])
        {
            wavelength_index = num_wavelength_points - 1;
        }
        else if (atmos.lambda[i] * 1e6 < wavelengths_in_microns[0])
        {
            wavelength_index = 0;
        }
        else
        {
            while (atmos.lambda[i] * 1e6 > wavelengths_in_microns[wavelength_index])
            {
                wavelength_index = wavelength_index+1;
            }
        }

        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                for(j=0; j<NTAU; j++)
                {
                    tau_em[l][m][j]=0.0;
                }
            }
        }

        /*Optical depth*/
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=90.0-PHASE && atmos.lon[m]<=270.0-PHASE)
                {
                    for(j=0; j<NTAU; j++)
                    {
                        Locate(NLAT, atmos.lat, theta_lat_solid[l][m][j], &o);
                        Locate(NLON, atmos.lon, phi_lon_solid[l][m][j]-PHASE, &c);

                        pressure = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.P_3d[o][c][j], atmos.P_3d[o][c+1][j], atmos.P_3d[o+1][c][j], atmos.P_3d[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);

                        //The pressure breaks this if it's too high
                        if (pressure > 9.99e9)
                        {
                            printf("Warning: pressures are too high\n");
                            pressure = 9.99e9;
                        }


                        if(atmos.T_3d[o][c][j] < 100.0 || atmos.T_3d[o][c+1][j] < 100.0)
                        {
                            temperature = 0.0;
                            temperature_3d[l][m][j] = temperature;
                        }

                        else
                        {
                            temperature = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.T_3d[o][c][j], atmos.T_3d[o][c+1][j], atmos.T_3d[o+1][c][j], atmos.T_3d[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            temperature_3d[l][m][j] = temperature;
                        }

                        Locate(NTEMP, opac.T, temperature, &g);
                        Locate(NPRESSURE, opac.P, pressure, &h);


                        /* Add doppler shift to signal, if turned on */
                        if(DOPPLER==1)
                        {
                            u_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ew[o][c][j], atmos.vel_ew[o][c+1][j], atmos.vel_ew[o+1][c][j], atmos.vel_ew[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ns[o][c][j], atmos.vel_ns[o][c+1][j], atmos.vel_ns[o+1][c][j], atmos.vel_ns[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            w_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ve[o][c][j], atmos.vel_ve[o][c+1][j], atmos.vel_ve[o+1][c][j], atmos.vel_ve[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_los = u_vel*sin(phi_lon_solid[l][m][j]*PI/180.0) + v_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*sin(theta_lat_solid[l][m][j]*PI/180.0) - w_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + cos(INPUT_INCLINATION)*omega*(R_PLANET + atmos.alt[j])*sin(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + R_VEL*cos((90.0-PHASE)*PI/180.0); /*Everything*/

                            delta_lam = atmos.lambda[i]*v_los/CLIGHT;
                            Locate(NLAMBDA, atmos.lambda, atmos.lambda[i]+delta_lam, &ii);

                            if(temperature < 100.0 || atmos.lambda[i]+delta_lam >= atmos.lambda[NLAMBDA-1] || atmos.lambda[i]+delta_lam <= atmos.lambda[0])
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint3D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1],
                                                  atmos.lambda[ii], atmos.lambda[ii+1],
                                                  opac.kappa[ii][h][g],
                                                  opac.kappa[ii][h][g+1],
                                                  opac.kappa[ii][h+1][g],
                                                  opac.kappa[ii][h+1][g+1],
                                                  opac.kappa[ii+1][h][g],
                                                  opac.kappa[ii+1][h][g+1],
                                                  opac.kappa[ii+1][h+1][g],
                                                  opac.kappa[ii+1][h+1][g+1],
                                                  temperature, pressure, atmos.lambda[i]+delta_lam);
                            }
                        }

                        /* Wind Only */
                        else if(DOPPLER==2)
                        {
                            u_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ew[o][c][j], atmos.vel_ew[o][c+1][j], atmos.vel_ew[o+1][c][j], atmos.vel_ew[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ns[o][c][j], atmos.vel_ns[o][c+1][j], atmos.vel_ns[o+1][c][j], atmos.vel_ns[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            w_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ve[o][c][j], atmos.vel_ve[o][c+1][j], atmos.vel_ve[o+1][c][j], atmos.vel_ve[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_los = u_vel*sin(phi_lon_solid[l][m][j]*PI/180.0) + v_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*sin(theta_lat_solid[l][m][j]*PI/180.0) - w_vel*cos(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + R_VEL*cos((90.0-PHASE)*PI/180.0);

                            delta_lam = atmos.lambda[i]*v_los/CLIGHT;
                            Locate(NLAMBDA, atmos.lambda, atmos.lambda[i]+delta_lam, &ii);

                            if(temperature < 250.0 || atmos.lambda[i]+delta_lam >= atmos.lambda[NLAMBDA-1] || atmos.lambda[i]+delta_lam <= atmos.lambda[0])
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint3D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1], atmos.lambda[ii], atmos.lambda[ii+1],
                                                  opac.kappa[ii][h][g],
                                                  opac.kappa[ii][h][g+1],
                                                  opac.kappa[ii][h+1][g],
                                                  opac.kappa[ii][h+1][g+1],
                                                  opac.kappa[ii+1][h][g],
                                                  opac.kappa[ii+1][h][g+1],
                                                  opac.kappa[ii+1][h+1][g],
                                                  opac.kappa[ii+1][h+1][g+1],
                                                  temperature, pressure, atmos.lambda[i]+delta_lam);
                            }
                        }

                        /* Rotation Only */
                        else if(DOPPLER==3){
                            u_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ew[o][c][j], atmos.vel_ew[o][c+1][j], atmos.vel_ew[o+1][c][j], atmos.vel_ew[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ns[o][c][j], atmos.vel_ns[o][c+1][j], atmos.vel_ns[o+1][c][j], atmos.vel_ns[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            w_vel = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], atmos.vel_ve[o][c][j], atmos.vel_ve[o][c+1][j], atmos.vel_ve[o+1][c][j], atmos.vel_ve[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            v_los = (cos(INPUT_INCLINATION)*(omega*(R_PLANET + atmos.alt[j])*sin(phi_lon_solid[l][m][j]*PI/180.0)*cos(theta_lat_solid[l][m][j]*PI/180.0) + R_VEL*cos((90.0-PHASE)*PI/180.0)));

                            delta_lam = atmos.lambda[i]*v_los/CLIGHT;
                            Locate(NLAMBDA, atmos.lambda, atmos.lambda[i]+delta_lam, &ii);

                            if(temperature < 250.0 || atmos.lambda[i]+delta_lam >= atmos.lambda[NLAMBDA-1] || atmos.lambda[i]+delta_lam <= atmos.lambda[0])
                            {
                                kappa_nu = 0.0;
                            }
                            else
                            {
                                kappa_nu = lint3D(opac.T[g], opac.T[g+1], opac.P[h],
                                                  opac.P[h+1], atmos.lambda[ii], atmos.lambda[ii+1],
                                                  opac.kappa[ii][h][g],
                                                  opac.kappa[ii][h][g+1],
                                                  opac.kappa[ii][h+1][g],
                                                  opac.kappa[ii][h+1][g+1],
                                                  opac.kappa[ii+1][h][g],
                                                  opac.kappa[ii+1][h][g+1],
                                                  opac.kappa[ii+1][h+1][g],
                                                  opac.kappa[ii+1][h+1][g+1],
                                                  temperature, pressure, atmos.lambda[i]+delta_lam);
                            }
                        }

                        /* No Doppler Effects at all */
                        else
                        {
                            if(temperature < 100.0)
                            {
                                kappa_nu = 0.0;
                            }

                            else
                            {
                                kappa_nu = lint2D(opac.T[g], opac.T[g+1],
                                                  opac.P[h], opac.P[h+1],
                                                  opac.kappa[i][h][g],
                                                  opac.kappa[i][h][g+1],
                                                  opac.kappa[i][h+1][g],
                                                  opac.kappa[i][h+1][g+1],
                                                  temperature, pressure);

                                //kappa_nu = lint2D(opac.T[g], opac.T[g+1],
                                //                  opac.P[h], opac.P[h+1],
                                //                  opac.kappa[100][h][g],
                                //                  opac.kappa[100][h][g+1],
                                //                  opac.kappa[100][h+1][g],
                                //                  opac.kappa[100][h+1][g+1],
                                //                  temperature, pressure);
                            }
                        }

                        kappa_nu_array[l][m][j] = kappa_nu;
                        dtau_em[l][m][j] = kappa_nu * dl[l][m][j];
                        pressure_array[l][m][j] = pressure;

                        if(CLOUDS==1)
                        {
                            // Find the nearest pressure index //
                            pressure_index = 0;
                            if (pressure > pressure_array_for_scattering_data_in_pascals[num_pressure_points-1])
                            {
                                pressure_index = num_pressure_points-1;
                            }
                            else
                            {
                                while (pressure > pressure_array_for_scattering_data_in_pascals[pressure_index])
                                {
                                    pressure_index = pressure_index+1;
                                }
                            }

                            aero_kappa_pre_qext_interp_1 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_1[o][c][j], aero_kappa_pre_qext_1[o][c+1][j], aero_kappa_pre_qext_1[o+1][c][j], aero_kappa_pre_qext_1[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_2 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_2[o][c][j], aero_kappa_pre_qext_2[o][c+1][j], aero_kappa_pre_qext_2[o+1][c][j], aero_kappa_pre_qext_2[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_3 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_3[o][c][j], aero_kappa_pre_qext_3[o][c+1][j], aero_kappa_pre_qext_3[o+1][c][j], aero_kappa_pre_qext_3[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_4 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_4[o][c][j], aero_kappa_pre_qext_4[o][c+1][j], aero_kappa_pre_qext_4[o+1][c][j], aero_kappa_pre_qext_4[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_5 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_5[o][c][j], aero_kappa_pre_qext_5[o][c+1][j], aero_kappa_pre_qext_5[o+1][c][j], aero_kappa_pre_qext_5[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_6 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_6[o][c][j], aero_kappa_pre_qext_6[o][c+1][j], aero_kappa_pre_qext_6[o+1][c][j], aero_kappa_pre_qext_6[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_7 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_7[o][c][j], aero_kappa_pre_qext_7[o][c+1][j], aero_kappa_pre_qext_7[o+1][c][j], aero_kappa_pre_qext_7[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_8 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_8[o][c][j], aero_kappa_pre_qext_8[o][c+1][j], aero_kappa_pre_qext_8[o+1][c][j], aero_kappa_pre_qext_8[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_9 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_9[o][c][j], aero_kappa_pre_qext_9[o][c+1][j], aero_kappa_pre_qext_9[o+1][c][j], aero_kappa_pre_qext_9[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_10 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_10[o][c][j], aero_kappa_pre_qext_10[o][c+1][j], aero_kappa_pre_qext_10[o+1][c][j], aero_kappa_pre_qext_10[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_11 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_11[o][c][j], aero_kappa_pre_qext_11[o][c+1][j], aero_kappa_pre_qext_11[o+1][c][j], aero_kappa_pre_qext_11[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_12 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_12[o][c][j], aero_kappa_pre_qext_12[o][c+1][j], aero_kappa_pre_qext_12[o+1][c][j], aero_kappa_pre_qext_12[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);
                            aero_kappa_pre_qext_interp_13 = lint2D(atmos.lon[c], atmos.lon[c+1], atmos.lat[o], atmos.lat[o+1], aero_kappa_pre_qext_13[o][c][j], aero_kappa_pre_qext_13[o][c+1][j], aero_kappa_pre_qext_13[o+1][c][j], aero_kappa_pre_qext_13[o+1][c+1][j], phi_lon_solid[l][m][j]-PHASE, theta_lat_solid[l][m][j]);

                            aero_kappa_1 = aero_kappa_pre_qext_interp_1 * KCl_wav_qext[pressure_index][0];
                            aero_kappa_2 = aero_kappa_pre_qext_interp_2 * ZnS_wav_qext[pressure_index][0];
                            aero_kappa_3 = aero_kappa_pre_qext_interp_3 * Na2S_wav_qext[pressure_index][0];
                            aero_kappa_4 = aero_kappa_pre_qext_interp_4 * MnS_wav_qext[pressure_index][0];
                            aero_kappa_5 = aero_kappa_pre_qext_interp_5 * Cr_wav_qext[pressure_index][0];
                            aero_kappa_6 = aero_kappa_pre_qext_interp_6 * SiO2_wav_qext[pressure_index][0];
                            aero_kappa_7 = aero_kappa_pre_qext_interp_7 * Mg2SiO4_wav_qext[pressure_index][0];
                            aero_kappa_8 = aero_kappa_pre_qext_interp_8 * VO_wav_qext[pressure_index][0];
                            aero_kappa_9 = aero_kappa_pre_qext_interp_9 * Ni_wav_qext[pressure_index][0];
                            aero_kappa_10 = aero_kappa_pre_qext_interp_10 * Fe_wav_qext[pressure_index][0];
                            aero_kappa_11 = aero_kappa_pre_qext_interp_11 * CaSiO4_wav_qext[pressure_index][0];
                            aero_kappa_12 = aero_kappa_pre_qext_interp_12 * CaTiO3_wav_qext[pressure_index][0];
                            aero_kappa_13 = aero_kappa_pre_qext_interp_13 * Al2O3_wav_qext[pressure_index][0];

                            // So all the cloud wavelength values are not wavelength dependant in the output files
                            // This takes the optical depth and adds the wavelength and particle size dependant scattering
                            // parameters. I think this might be good code except it should all be arrays
                            total_cloud_kappa = aero_kappa_1 + aero_kappa_2 + aero_kappa_3 + aero_kappa_4 + \
                                                aero_kappa_5 + aero_kappa_6 + aero_kappa_7 + aero_kappa_8 + \
                                                aero_kappa_9 + aero_kappa_10 + aero_kappa_11 + aero_kappa_12 + \
                                                aero_kappa_13;

                            if (dtau_em[l][m][j] < 1e-3 || total_cloud_kappa < 1e-5 || kappa_nu < 1e-10 || dl[l][m][j] < 1e-10)
                            //if (dtau_em[l][m][j] < 1e-50 || total_cloud_kappa < 1e-50 || kappa_nu < 1e-50 || dl[l][m][j] < 1e-50)
                            {
                                pi0_tot[l][m][j] = 0.0;
                                asym_tot[l][m][j] = 0.0;
                                dtau_em[l][m][j] = kappa_nu * dl[l][m][j];
                            }
                            else
                            {
                                dtau_em[l][m][j] = (kappa_nu_array[l][m][j] + total_cloud_kappa) * dl[l][m][j];
                                temp_value = dl[l][m][j] / dtau_em[l][m][j];

                                weight_1 = aero_kappa_1 * temp_value;
                                weight_2 = aero_kappa_2 * temp_value;
                                weight_3 = aero_kappa_3 * temp_value;
                                weight_4 = aero_kappa_4 * temp_value;
                                weight_5 = aero_kappa_5 * temp_value;
                                weight_6 = aero_kappa_6 * temp_value;
                                weight_7 = aero_kappa_7 * temp_value;
                                weight_8 = aero_kappa_8 * temp_value;
                                weight_9 = aero_kappa_9 * temp_value;
                                weight_10 = aero_kappa_10 * temp_value;
                                weight_11 = aero_kappa_11 * temp_value;
                                weight_12 = aero_kappa_12 * temp_value;
                                weight_13 = aero_kappa_13 * temp_value;

                                pi0_tot[l][m][j] =  (weight_1  * KCl_wav_pi0[pressure_index][0] + \
                                                     weight_2  * ZnS_wav_pi0[pressure_index][0] + \
                                                     weight_3  * Na2S_wav_pi0[pressure_index][0] + \
                                                     weight_4  * MnS_wav_pi0[pressure_index][0] + \
                                                     weight_5  * Cr_wav_pi0[pressure_index][0] + \
                                                     weight_6  * SiO2_wav_pi0[pressure_index][0] + \
                                                     weight_7  * Mg2SiO4_wav_pi0[pressure_index][0] + \
                                                     weight_8  * VO_wav_pi0[pressure_index][0] + \
                                                     weight_9  * Ni_wav_pi0[pressure_index][0] + \
                                                     weight_10 * Fe_wav_pi0[pressure_index][0] + \
                                                     weight_11 * CaSiO4_wav_pi0[pressure_index][0] + \
                                                     weight_12 * CaTiO3_wav_pi0[pressure_index][0] + \
                                                     weight_13 * Al2O3_wav_pi0[pressure_index][0]);

                                asym_tot[l][m][j] = (weight_1  * KCl_wav_gg[pressure_index][0] + \
                                                     weight_2  * ZnS_wav_gg[pressure_index][0]  + \
                                                     weight_3  * Na2S_wav_gg[pressure_index][0] + \
                                                     weight_4  * MnS_wav_gg[pressure_index][0] + \
                                                     weight_5  * Cr_wav_gg[pressure_index][0] + \
                                                     weight_6  * SiO2_wav_gg[pressure_index][0] + \
                                                     weight_7  * Mg2SiO4_wav_gg[pressure_index][0] + \
                                                     weight_8  * VO_wav_gg[pressure_index][0] + \
                                                     weight_9  * Ni_wav_gg[pressure_index][0] + \
                                                     weight_10 * Fe_wav_gg[pressure_index][0] + \
                                                     weight_11 * CaSiO4_wav_gg[pressure_index][0] + \
                                                     weight_12 * CaTiO3_wav_gg[pressure_index][0] + \
                                                     weight_13 * Al2O3_wav_gg[pressure_index][0]);
                            }
                        }
                        // if clouds are turned off, need to set scattering params to zero
                        else
                        {
                            pi0_tot[l][m][j] = 0.0;
                            asym_tot[l][m][j] = 0.0;
                        }

                    }
                }
            }
        }

        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=90.0-PHASE && atmos.lon[m]<=270.0-PHASE)
                {
                    for(j=0; j<NTAU; j++)
                    {
                        if(j==0)
                        {
                            tau_em[l][m][j]=dtau_em[l][m][j];
                        }
                        else if (j!=0)
                        {
                            tau_em[l][m][j] = tau_em[l][m][j-1] + dtau_em[l][m][j];
                        }
                    }
                }
            }
        }


        //Calculate the intensity of emergent rays at each latitude and longitude

        // C HARADA -- update for 2stream //
        running_sum = 0.0;
        average = 0.0;
        good_intensity = 0.0;

        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                intensity[l][m] = 0;
                if(atmos.lon[m]>=90.0-PHASE && atmos.lon[m]<=270.0-PHASE)
                {

                    // Find min vert level for 2stream //
                    kmin = 0;
                    for (j = 0; j<NTAU; j++)
                    {
                        if (dtau_em[l][m][j] < 1e-10 || tau_em[l][m][j] < 1e-10 || temperature_3d[l][m][j] < 250)
                        {
                          kmin = j+1;
                        }
                    }

                    if (kmin >= NTAU)
                    {
                        intensity[l][m] = 0;
                    }

                    if (atmos.incident_frac[l][m][NTAU-10] < 0)
                    {
                        atmos.incident_frac[l][m][NTAU-10] = 0;
                    }

                    else
                    {
                        intensity[l][m] = two_stream(NTAU, kmin, pi0_tot[l][m], asym_tot[l][m], temperature_3d[l][m], tau_em[l][m], \
                                                     CLIGHT / atmos.lambda[i], CLIGHT / atmos.lambda[i] - CLIGHT / atmos.lambda[i+1],
                                                     atmos.incident_frac[l][m][NTAU-10], dtau_em[l][m]);
                    }
                }
            }
        }




        /*
        // ~~~ THIS IS THE OLD RT ROUTINE ~~~ //
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=90.0-PHASE && atmos.lon[m]<=270.0-PHASE)
                {
                    intensity[l][m] = Planck(atmos.T_3d[l][m][NTAU-1], atmos.lambda[i]) * exp(-tau_em[l][m][NTAU-1]);
                }
            }
        }

        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=90.0-PHASE && atmos.lon[m]<=270.0-PHASE)
                {
                    for(j=0; j<NTAU; j++)
                    {

                        intensity[l][m] += Planck(temperature_3d[l][m][j], atmos.lambda[i]) * exp(-tau_em[l][m][j]) * dtau_em[l][m][j];
                    }
                }
            }
        }
        */





        /*Calculate the total flux received by us*/
        flux_pl[i] = 0.0;
        for(l=0; l<NLAT; l++)
        {
            for(m=0; m<NLON; m++)
            {
                if(atmos.lon[m]>=90.0-PHASE && atmos.lon[m]<=270.0-PHASE)
                {
                    flux_pl[i] += intensity[l][m] * SQ(cos(theta[l][m])) * cos(phi[l][m]-PI) * dtheta[l][m] * dphi[l][m];
                }
            }
        }
        //fclose(fptr);


        if(i % 100 == 0)
        {
            printf("%d out of %d lines (phase: %06.2f)\n", i, NLAMBDA, PHASE);
        }

        fprintf(file, "%10.8le\t%le\n", atmos.lambda[i], flux_pl[i] * PI/solid);

    }

    fclose(file);

    return 0;
}