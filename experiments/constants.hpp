#define MINUS_ONE (-1.0)
#define PLUS_ONE 1.0
#define BETA_LO 0.9401653
#define BETA_HI 0.9401658
#define LOWER_CUTOFF_RATIO 0.94016566
#define UPPER_CUTOFF_RATIO 0.9401658
#define STEP_1_RAD 0.0099999999
#define EPS_LIMIT 0.0000001
#define HARD_LO_B1 (-0.16247834)
#define HARD_LO_B2 HARD_LO_B1
#define HARD_LO_B12 (-0.67504335)
#define HARD_HI_B1 0.16247834
#define HARD_HI_B2 HARD_HI_B1
#define HARD_HI_B12 (-0.67504335)
#define HARD_B12 HARD_LO_B12
#define STEP_3_RAD 0.05
#define T_MAX 0.1

#define TYPE_3_EPS 0.000001
#define TYPE_3_B1_HARD (-0.1824167935)
#define TYPE_3_B2_HARD (-0.1824167935)
#define TYPE_3_BETA 0.9539799006
#define TYPE_3_COARSE_BETA_LO 0.95
#define TYPE_3_COARSE_BETA_HI 0.96
#define TYPE_3_FINE_BETA_LO 0.9539798
#define TYPE_3_FINE_BETA_HI 0.95398
#define TYPE_3_LOWER_BOUND 0.001

#define ALPHA_HI 0.95
#define ALPHA_LO 0.9461
#define OBJ_HI ((ALPHA_HI - 1) / ALPHA_HI)
#define OBJ_LO ((ALPHA_LO - 1) / ALPHA_LO)
#define TYPE_4_EPS EPS_LIMIT
#define TYPE_4_HARD_EPS 0.000001
#define TYPE_4_HARD_EPS_ALT 0.01
#define TYPE_4_HARD_EPS_ALT2 0.0001
#define TYPE_4_HARD_BOUND 0.8284084
#define TYPE_4_B1_HARD 0.1489442419
#define TYPE_4_B2_HARD 0.1489442419
#define TYPE_4_T1 ((1 - TYPE_4_B1_HARD) / 2)
#define TYPE_4_T2 ((1 + TYPE_4_B1_HARD) / 2)
#define TYPE_4_T_EPS 0.01
#define TYPE_4_P1 (Arb(0.085795671264819))
#define TYPE_4_P2 (Arb(0.085795671264819))
#define TYPE_4_P3 (Arb(0.173738866098353))
#define TYPE_4_P4 (Arb(0.483078448842371))
#define TYPE_4_P5 (Arb(0.085795671264819))
#define TYPE_4_P6 (Arb(0.085795671264819))
#define TYPE_4_RATIO 0.9461598105

#define TYPE_5_COARSE_BETA_LO .94
#define TYPE_5_COARSE_BETA_HI .9405
#define TYPE_5_BETA_LO 0.9401653
#define TYPE_5_BETA_HI 0.9401658
#define TYPE_5_B1_HARD 0.16247834
#define TYPE_5_B2_HARD 0.16247834
#define TYPE_5_EPS 0.000001
#define TYPE_5_LOWER_BOUND 0.001
//#define TYPE_5_RATIO_UPPER 0.941

//#define DEBUG
