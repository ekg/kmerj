#include "melting.hpp"

namespace kmerj {

namespace melting {

// adapted from https://github.com/Sely85/dna_melting

double wallace_rule(double sequence_length, double a_count , double c_count, double g_count, double t_count)  {
    //Marmur and Doty, 1962
    double wallace_melting_temperature;
    if (sequence_length <= 15) {
        wallace_melting_temperature=2*(a_count+t_count)+4*(c_count+g_count);
    } else {
        wallace_melting_temperature=69.3+(((41*(c_count+g_count))/sequence_length)-(650/sequence_length));
    }
    return wallace_melting_temperature;
}

double salt(double salt_conc, double a_count , double c_count, double g_count, double t_count) {
    //Howley et al., 1979

    double salt_adjusted_melting_temperature;

    salt_adjusted_melting_temperature=81.5+16.6*log10(salt_conc)+41.0*((c_count+g_count)/(a_count+c_count+g_count+t_count))-(675.0/(a_count+c_count+g_count+t_count));
    
    return salt_adjusted_melting_temperature;
}

double khandelwal(int sequence_length, const std::string& specie, double salt_conc, double dna_conc) {
    //Khandelwal and Bhyravabhotla, 2010

    double khandelwal_melting_temperature;

    double ee = 0;
    double strength = 0;
    int i = 0;

    double gc=13;
    double ac=10;
    double gt=10;
    double at=7;
    double cc=11;
    double tc=8;
    double ct=8;
    double tt=5;
    double gg=11;
    double ag=8;
    double ga=8;
    double aa=5;
    double cg=10;
    double tg=7;
    double ca=7;
    double ta=4;

    i=0;
    ee=0;

    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') strength=strength+aa;
        if (specie[i]=='A' && specie[i+1]=='C') strength=strength+ac;
        if (specie[i]=='A' && specie[i+1]=='G') strength=strength+ag;
        if (specie[i]=='A' && specie[i+1]=='T') strength=strength+at;

        if (specie[i]=='C' && specie[i+1]=='A') strength=strength+ca;
        if (specie[i]=='C' && specie[i+1]=='C') strength=strength+cc;
        if (specie[i]=='C' && specie[i+1]=='G') strength=strength+cg;
        if (specie[i]=='C' && specie[i+1]=='T') strength=strength+ct;

        if (specie[i]=='G' && specie[i+1]=='A') strength=strength+ga;
        if (specie[i]=='G' && specie[i+1]=='C') strength=strength+gc;
        if (specie[i]=='G' && specie[i+1]=='G') strength=strength+gg;
        if (specie[i]=='G' && specie[i+1]=='T') strength=strength+gt;

        if (specie[i]=='T' && specie[i+1]=='A') strength=strength+ta;
        if (specie[i]=='T' && specie[i+1]=='C') strength=strength+tc;
        if (specie[i]=='T' && specie[i+1]=='G') strength=strength+tg;
        if (specie[i]=='T' && specie[i+1]=='T') strength=strength+tt;
    }

    ee = strength/sequence_length;

    khandelwal_melting_temperature=7.35*ee+17.34*log(sequence_length)+4.96*log(salt_conc)+0.89*log(dna_conc)-25.42;

    return khandelwal_melting_temperature;

}

double bre_nearest_neighbor(int sequence_length, const std::string& specie, double salt_conc, double dna_conc, double a_count, double c_count, double g_count, double t_count)
{
    //Breslauer, Frank, Blocker and Marky, 1986

    double nearest_neighbor_melting_temperature;

    int i;
    double simm_corr;
    double non_self_compl;
    double only_at;
    double any_cg;          

    //deltaS
    simm_corr=-1.34;	 //kcal/(K mol)
    non_self_compl=0.0;  //kcal/(K mol)
    only_at=-20.13;	 //kcal/(K mol)
    any_cg=-16.77;          //kcal/(K mol)

    double h_aa, h_ac, h_ag, h_at, h_ca, h_cc, h_cg, h_ct, h_ga, h_gc, h_gg, h_gt, h_ta, h_tc, h_tg, h_tt;
    double s_aa, s_ac, s_ag, s_at, s_ca, s_cc, s_cg, s_ct, s_ga, s_gc, s_gg, s_gt, s_ta, s_tc, s_tg, s_tt;
    double deltah_d, deltas_d, deltas_self;

    deltah_d=0;
    deltas_d=0;

    i=0;

    //Enthaply
    //kcal/mol
    h_aa=h_tt=-9.1;
    h_at=-8.6;
    h_ta=-6.0;
    h_ca=h_tg=-5.8;
    h_gt=h_ac=-6.5;
    h_ct=h_ag=-7.8;
    h_ga=h_tc=-5.6;
    h_cg=-11.9;
    h_gc=-11.1;
    h_gg=h_cc=-11.0;


    //Entropy
    //cal/(K mol)
    s_aa=s_tt=-24.0;
    s_at=-23.9;
    s_ta=-16.9;
    s_ca=s_tg=-12.9;
    s_gt=s_ac=-17.3;
    s_ct=s_ag=-20.8;
    s_ga=s_tc=-13.5;
    s_cg=-27.8;
    s_gc=-26.7;
    s_gg=s_cc=-26.6;


//ENTHALPY
//kcal/mol
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltah_d=deltah_d+h_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltah_d=deltah_d+h_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltah_d=deltah_d+h_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltah_d=deltah_d+h_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltah_d=deltah_d+h_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltah_d=deltah_d+h_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltah_d=deltah_d+h_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltah_d=deltah_d+h_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltah_d=deltah_d+h_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltah_d=deltah_d+h_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltah_d=deltah_d+h_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltah_d=deltah_d+h_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltah_d=deltah_d+h_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltah_d=deltah_d+h_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltah_d=deltah_d+h_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltah_d=deltah_d+h_tt;
    }


    //ENTROPY
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltas_d=deltas_d+s_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltas_d=deltas_d+s_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltas_d=deltas_d+s_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltas_d=deltas_d+s_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltas_d=deltas_d+s_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltas_d=deltas_d+s_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltas_d=deltas_d+s_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltas_d=deltas_d+s_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltas_d=deltas_d+s_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltas_d=deltas_d+s_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltas_d=deltas_d+s_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltas_d=deltas_d+s_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltas_d=deltas_d+s_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltas_d=deltas_d+s_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltas_d=deltas_d+s_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltas_d=deltas_d+s_tt;

    }

    double R=1.987; //cal/(K mol)

    double b;
    if (specie[i]==specie[sequence_length] && specie[i+1]==specie[sequence_length-1] && specie[i+2]==specie[sequence_length-2]  && specie[i+3]==specie[sequence_length-3] && specie[i+4]==specie[sequence_length-4] && specie[i+5]==specie[sequence_length-5]) b=1;
    else b=4;

    if (specie[i]==specie[sequence_length] && specie[i+1]==specie[sequence_length-1] && specie[i+2]==specie[sequence_length-2]  && specie[i+3]==specie[sequence_length-3] && specie[i+4]==specie[sequence_length-4] && specie[i+5]==specie[sequence_length-5]) deltas_self=simm_corr;


    else deltas_self=non_self_compl;

    double deltah_i=0;
    double deltas_i;

    if (c_count==0 && g_count==0) deltas_i=only_at;
    else deltas_i=any_cg;

    //enthalpy&R in cal 
    double num=deltah_d*1000+deltah_i*1000;
    double den=deltas_d+deltas_i+deltas_self+R*log(dna_conc/b);
    double salt_adj=16.6*log10(salt_conc);

    nearest_neighbor_melting_temperature=num/den+salt_adj;

    return nearest_neighbor_melting_temperature;

}



double san_nearest_neighbor(int sequence_length, const std::string& specie, double salt_conc, double dna_conc, double a_count, double c_count, double g_count, double t_count)
{
    //SantaLucia, Allawi and Seneviratne, 1996

    double nearest_neighbor_melting_temperature;

    int i;
    double simm_corr;
    double non_self_compl;
    double only_at;
    double any_cg;          



    //deltaS
    simm_corr=-1.4;	 //kcal/(K mol)
    non_self_compl=0.0;  //kcal/(K mol)
    only_at=-9.0;	 //kcal/(K mol)
    any_cg=-5.9;          //kcal/(K mol)


    double h_aa, h_ac, h_ag, h_at, h_ca, h_cc, h_cg, h_ct, h_ga, h_gc, h_gg, h_gt, h_ta, h_tc, h_tg, h_tt;
    double s_aa, s_ac, s_ag, s_at, s_ca, s_cc, s_cg, s_ct, s_ga, s_gc, s_gg, s_gt, s_ta, s_tc, s_tg, s_tt;
    double deltah_d, deltas_d, deltas_self;

    deltah_d=0;
    deltas_d=0;

    i=0;

    //Enthalpy
    //kcal/(mol
    h_aa=h_tt=-8.4;
    h_at=-6.5;
    h_ta=-6.3;
    h_ca=h_tg=-7.4;
    h_gt=h_ac=-8.6;
    h_ct=h_ag=-6.1;
    h_ga=h_tc=-7.7;
    h_cg=-10.1;
    h_gc=-11.1;
    h_gg=h_cc=-6.7;


    //Entropy
    //cal/(K mol)
    s_aa=s_tt=-23.6;
    s_at=-18.8;
    s_ta=-18.5;
    s_ca=s_tg=-19.3;
    s_gt=s_ac=-23.0;
    s_ct=s_ag=-16.1;
    s_ga=s_tc=-20.3;
    s_cg=-25.5;
    s_gc=-28.4;
    s_gg=s_cc=-15.6;


//ENTHALPY
//kcal/mol
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltah_d=deltah_d+h_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltah_d=deltah_d+h_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltah_d=deltah_d+h_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltah_d=deltah_d+h_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltah_d=deltah_d+h_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltah_d=deltah_d+h_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltah_d=deltah_d+h_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltah_d=deltah_d+h_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltah_d=deltah_d+h_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltah_d=deltah_d+h_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltah_d=deltah_d+h_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltah_d=deltah_d+h_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltah_d=deltah_d+h_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltah_d=deltah_d+h_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltah_d=deltah_d+h_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltah_d=deltah_d+h_tt;
    }


    //ENTROPY
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltas_d=deltas_d+s_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltas_d=deltas_d+s_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltas_d=deltas_d+s_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltas_d=deltas_d+s_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltas_d=deltas_d+s_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltas_d=deltas_d+s_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltas_d=deltas_d+s_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltas_d=deltas_d+s_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltas_d=deltas_d+s_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltas_d=deltas_d+s_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltas_d=deltas_d+s_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltas_d=deltas_d+s_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltas_d=deltas_d+s_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltas_d=deltas_d+s_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltas_d=deltas_d+s_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltas_d=deltas_d+s_tt;

    }

    double R=1.987; //cal/(K mol)

    double b;
    if (specie[i]==specie[sequence_length] && specie[i+1]==specie[sequence_length-1] && specie[i+2]==specie[sequence_length-2]  && specie[i+3]==specie[sequence_length-3] && specie[i+4]==specie[sequence_length-4] && specie[i+5]==specie[sequence_length-5]) b=1;
    else b=4;

    if (specie[i]==specie[sequence_length] && specie[i+1]==specie[sequence_length-1] && specie[i+2]==specie[sequence_length-2]  && specie[i+3]==specie[sequence_length-3] && specie[i+4]==specie[sequence_length-4] && specie[i+5]==specie[sequence_length-5]) deltas_self=simm_corr;


    else deltas_self=non_self_compl;

    double deltah_i=0;
    double deltas_i;

    if (c_count==0 && g_count==0) deltas_i=only_at;
    else deltas_i=any_cg;

    //enthalpy&R in cal 
    double num=deltah_d*1000+deltah_i*1000;
    double den=deltas_d+deltas_i+deltas_self+R*log(dna_conc/b);
    double salt_adj=16.6*log10(salt_conc);

    nearest_neighbor_melting_temperature=num/den+salt_adj;

    return nearest_neighbor_melting_temperature;

}




double sug_nearest_neighbor(int sequence_length, const std::string& specie, double salt_conc, double dna_conc, double a_count, double c_count, double g_count, double t_count)
{
    //Sugimoto, Nakano, Yoneyama and Honda, 1996

    double nearest_neighbor_melting_temperature;

    int i;
    double simm_corr;
    double non_self_compl;
    double only_at;
    double any_cg;          


    //deltaS
    simm_corr=-1.4;	 //kcal/(K mol)
    non_self_compl=0.0;  //kcal/(K mol)
    only_at=-9.0;	 //kcal/(K mol)
    any_cg=-9.0;          //kcal/(K mol)


    double h_aa, h_ac, h_ag, h_at, h_ca, h_cc, h_cg, h_ct, h_ga, h_gc, h_gg, h_gt, h_ta, h_tc, h_tg, h_tt;
    double s_aa, s_ac, s_ag, s_at, s_ca, s_cc, s_cg, s_ct, s_ga, s_gc, s_gg, s_gt, s_ta, s_tc, s_tg, s_tt;
    double deltah_d, deltas_d, deltas_self;

    deltah_d=0;
    deltas_d=0;

    i=0;

    //Enthalpy
    //kcal/(mol
    h_aa=h_tt=-8.0;
    h_at=-5.6;
    h_ta=-6.6;
    h_ca=h_tg=-8.2;
    h_gt=h_ac=-9.4;
    h_ct=h_ag=-6.6;
    h_ga=h_tc=-8.8;
    h_cg=-11.8;
    h_gc=-10.5;
    h_gg=h_cc=-10.9;


    //Entropy
    //cal/(K mol)
    s_aa=s_tt=-21.9;
    s_at=-15.2;
    s_ta=-18.4;
    s_ca=s_tg=-21.0;
    s_gt=s_ac=-25.5;
    s_ct=s_ag=-16.4;
    s_ga=s_tc=-23.5;
    s_cg=-29.0;
    s_gc=-26.4;
    s_gg=s_cc=-28.4;


//ENTHALPY
//kcal/mol
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltah_d=deltah_d+h_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltah_d=deltah_d+h_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltah_d=deltah_d+h_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltah_d=deltah_d+h_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltah_d=deltah_d+h_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltah_d=deltah_d+h_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltah_d=deltah_d+h_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltah_d=deltah_d+h_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltah_d=deltah_d+h_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltah_d=deltah_d+h_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltah_d=deltah_d+h_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltah_d=deltah_d+h_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltah_d=deltah_d+h_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltah_d=deltah_d+h_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltah_d=deltah_d+h_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltah_d=deltah_d+h_tt;
    }

    //ENTROPY
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltas_d=deltas_d+s_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltas_d=deltas_d+s_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltas_d=deltas_d+s_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltas_d=deltas_d+s_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltas_d=deltas_d+s_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltas_d=deltas_d+s_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltas_d=deltas_d+s_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltas_d=deltas_d+s_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltas_d=deltas_d+s_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltas_d=deltas_d+s_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltas_d=deltas_d+s_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltas_d=deltas_d+s_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltas_d=deltas_d+s_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltas_d=deltas_d+s_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltas_d=deltas_d+s_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltas_d=deltas_d+s_tt;

    }

    double R=1.987; //cal/(K mol)

    double b;
    if (specie[i]==specie[sequence_length] && specie[i+1]==specie[sequence_length-1] && specie[i+2]==specie[sequence_length-2]  && specie[i+3]==specie[sequence_length-3] && specie[i+4]==specie[sequence_length-4] && specie[i+5]==specie[sequence_length-5]) b=1;
    else b=4;

    if (specie[i]==specie[sequence_length] && specie[i+1]==specie[sequence_length-1] && specie[i+2]==specie[sequence_length-2]  && specie[i+3]==specie[sequence_length-3] && specie[i+4]==specie[sequence_length-4] && specie[i+5]==specie[sequence_length-5]) deltas_self=simm_corr;


    else deltas_self=non_self_compl;

    double deltah_i=0;
    double deltas_i;

    if (c_count==0 && g_count==0) deltas_i=only_at;
    else deltas_i=any_cg;

    //enthalpy&R in cal 
    double num=deltah_d*1000+deltah_i*1000;
    double den=deltas_d+deltas_i+deltas_self+R*log(dna_conc/b);
    double salt_adj=16.6*log10(salt_conc);

    nearest_neighbor_melting_temperature=num/den+salt_adj;

    return nearest_neighbor_melting_temperature;

}




double consensus(int sequence_length, const std::string& specie, double salt_conc, double dna_conc, double a_count, double c_count, double g_count, double t_count)
{
    //Panjkovich and Melo, 2005

    double consensus_melting_temperature;

    int i;

    double bre_tm = bre_nearest_neighbor(sequence_length, specie, salt_conc, dna_conc, a_count, c_count, g_count, t_count);
    double san_tm = san_nearest_neighbor(sequence_length, specie, salt_conc, dna_conc, a_count, c_count, g_count, t_count);
    double sug_tm = sug_nearest_neighbor(sequence_length, specie, salt_conc, dna_conc, a_count, c_count, g_count, t_count);


    double gc_content = (c_count + g_count)/(a_count + c_count + g_count + t_count)*100;


    if (sequence_length>=16 && sequence_length<=18 && gc_content>=30 && gc_content<=50){
        consensus_melting_temperature=(bre_tm+san_tm+sug_tm)/3;
        std::cout << "Full consensus sequence " << std::endl;
    }

    if (sequence_length>=21 && sequence_length<=22 && gc_content>=0 && gc_content<=10 || sequence_length>=16 && sequence_length<=29 && gc_content>=10 && gc_content<=20 || sequence_length>=16 && sequence_length<=26 && gc_content>=20 && gc_content<=30 || sequence_length>=19 && sequence_length<=21 && gc_content>=30 && gc_content<=40){
        consensus_melting_temperature=(bre_tm+sug_tm)/2;
        std::cout << "Bre&Sug consensus sequence " << std::endl;
    }

    if (sequence_length>=19 && sequence_length<=30 && gc_content>=40 && gc_content<=50 || sequence_length>=16 && sequence_length<=30 && gc_content>=50 && gc_content<=80 || sequence_length==16 && gc_content>=80 && gc_content<=90 || sequence_length==17 && gc_content>=80 && gc_content<=90 || sequence_length==20 && gc_content>=80 && gc_content<=90){
        consensus_melting_temperature=(san_tm+sug_tm)/2;
        std::cout << "San&Sug consensus sequence " << std::endl;
    }
    
    else 
    {
        consensus_melting_temperature=(bre_tm+san_tm+sug_tm)/3;
        std::cout << "Non-consensus sequence " << std::endl;
    }

    return consensus_melting_temperature;

}


double bre_enthalpy(int sequence_length, const std::string& specie)
{
    double bre_enthalpy;
    double deltah;
    int i;

    deltah=0;
    i=0;

    double h_aa, h_ac, h_ag, h_at, h_ca, h_cc, h_cg, h_ct, h_ga, h_gc, h_gg, h_gt, h_ta, h_tc, h_tg, h_tt;



    //Enthaply
    //kcal/mol
    h_aa=h_tt=-9.1;
    h_at=-8.6;
    h_ta=-6.0;
    h_ca=h_tg=-5.8;
    h_gt=h_ac=-6.5;
    h_ct=h_ag=-7.8;
    h_ga=h_tc=-5.6;
    h_cg=-11.9;
    h_gc=-11.1;
    h_gg=h_cc=-11.0;



//ENTHALPY
//kcal/mol
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltah=deltah+h_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltah=deltah+h_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltah=deltah+h_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltah=deltah+h_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltah=deltah+h_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltah=deltah+h_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltah=deltah+h_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltah=deltah+h_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltah=deltah+h_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltah=deltah+h_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltah=deltah+h_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltah=deltah+h_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltah=deltah+h_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltah=deltah+h_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltah=deltah+h_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltah=deltah+h_tt;
    }


    deltah=deltah*1000;
    bre_enthalpy=deltah;

    return bre_enthalpy;

}



double bre_entropy(int sequence_length, const std::string& specie)
{

    double bre_entropy;

    int i;
    double s_aa, s_ac, s_ag, s_at, s_ca, s_cc, s_cg, s_ct, s_ga, s_gc, s_gg, s_gt, s_ta, s_tc, s_tg, s_tt;
    double deltas;

    deltas=0;
    i=0;

    //Entropy
    //cal/(K mol)
    s_aa=s_tt=-24.0;
    s_at=-23.9;
    s_ta=-16.9;
    s_ca=s_tg=-12.9;
    s_gt=s_ac=-17.3;
    s_ct=s_ag=-20.8;
    s_ga=s_tc=-13.5;
    s_cg=-27.8;
    s_gc=-26.7;
    s_gg=s_cc=-26.6;


    //ENTROPY
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltas=deltas+s_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltas=deltas+s_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltas=deltas+s_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltas=deltas+s_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltas=deltas+s_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltas=deltas+s_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltas=deltas+s_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltas=deltas+s_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltas=deltas+s_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltas=deltas+s_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltas=deltas+s_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltas=deltas+s_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltas=deltas+s_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltas=deltas+s_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltas=deltas+s_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltas=deltas+s_tt;

    }

    bre_entropy=deltas;

    return bre_entropy;

}




double bre_melting_curve(int sequence_length, const std::string& specie, double dna_conc)
{
    double f, keq, ctkeq;

    double deltah = bre_enthalpy(sequence_length, specie);
    double deltas = bre_entropy(sequence_length, specie);

    double t;
    double R=1.987; //cal/(K mol)

    t=0.0;
    std::ofstream fileout("bre_melting_curve.out");


    for (t=0.0; t<=700.0; t=t+0.5){
        ctkeq=dna_conc*exp((deltas/R)-(deltah/(R*t)));
        f=(1+ctkeq-sqrt(1+2*ctkeq))/ctkeq;
        fileout << t << " " << f << std::endl;    
    }

    return 0;

}




double san_enthalpy(int sequence_length, const std::string& specie)
{
    double san_enthalpy;
    int i;

    double h_aa, h_ac, h_ag, h_at, h_ca, h_cc, h_cg, h_ct, h_ga, h_gc, h_gg, h_gt, h_ta, h_tc, h_tg, h_tt;
    double deltah;

    deltah=0;
    i=0;

    //Enthalpy
    //kcal/(mol
    h_aa=h_tt=-8.4;
    h_at=-6.5;
    h_ta=-6.3;
    h_ca=h_tg=-7.4;
    h_gt=h_ac=-8.6;
    h_ct=h_ag=-6.1;
    h_ga=h_tc=-7.7;
    h_cg=-10.1;
    h_gc=-11.1;
    h_gg=h_cc=-6.7;


//ENTHALPY
//kcal/mol
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltah=deltah+h_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltah=deltah+h_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltah=deltah+h_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltah=deltah+h_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltah=deltah+h_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltah=deltah+h_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltah=deltah+h_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltah=deltah+h_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltah=deltah+h_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltah=deltah+h_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltah=deltah+h_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltah=deltah+h_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltah=deltah+h_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltah=deltah+h_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltah=deltah+h_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltah=deltah+h_tt;
    }

    san_enthalpy=deltah*1000;
    return san_enthalpy;

}

double san_entropy(int sequence_length, const std::string& specie)
{
    double san_entropy;
    int i;

    double s_aa, s_ac, s_ag, s_at, s_ca, s_cc, s_cg, s_ct, s_ga, s_gc, s_gg, s_gt, s_ta, s_tc, s_tg, s_tt;
    double deltas;

    deltas=0;
    i=0;


    //Entropy
    //cal/(K mol)
    s_aa=s_tt=-23.6;
    s_at=-18.8;
    s_ta=-18.5;
    s_ca=s_tg=-19.3;
    s_gt=s_ac=-23.0;
    s_ct=s_ag=-16.1;
    s_ga=s_tc=-20.3;
    s_cg=-25.5;
    s_gc=-28.4;
    s_gg=s_cc=-15.6;


    //ENTROPY
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltas=deltas+s_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltas=deltas+s_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltas=deltas+s_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltas=deltas+s_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltas=deltas+s_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltas=deltas+s_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltas=deltas+s_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltas=deltas+s_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltas=deltas+s_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltas=deltas+s_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltas=deltas+s_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltas=deltas+s_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltas=deltas+s_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltas=deltas+s_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltas=deltas+s_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltas=deltas+s_tt;

    }

    san_entropy=deltas;
    return san_entropy;

}





double san_melting_curve(int sequence_length, const std::string& specie, double dna_conc)
{
    double f, keq, ctkeq;

    double deltah = san_enthalpy(sequence_length, specie);
    double deltas = san_entropy(sequence_length, specie);

    double t;
    double R=1.987; //cal/(K mol)

    t=0.0;
    std::ofstream fileout("san_melting_curve.out");


    for (t=0.0; t<=700.0; t=t+0.5){
        ctkeq=dna_conc*exp((deltas/R)-(deltah/(R*t)));
        f=(1+ctkeq-sqrt(1+2*ctkeq))/ctkeq;
        fileout << t << " " << f << std::endl;    
    }

    return 0;

}




double sug_enthalpy(int sequence_length, const std::string& specie)
{
    double sug_enthalpy;
    int i;
    double h_aa, h_ac, h_ag, h_at, h_ca, h_cc, h_cg, h_ct, h_ga, h_gc, h_gg, h_gt, h_ta, h_tc, h_tg, h_tt;
    double deltah;

    deltah=0;
    i=0;

    //Enthalpy
    //kcal/mol
    h_aa=h_tt=-8.0;
    h_at=-5.6;
    h_ta=-6.6;
    h_ca=h_tg=-8.2;
    h_gt=h_ac=-9.4;
    h_ct=h_ag=-6.6;
    h_ga=h_tc=-8.8;
    h_cg=-11.8;
    h_gc=-10.5;
    h_gg=h_cc=-10.9;

    //ENTHALPY
    //kcal/mol
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltah=deltah+h_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltah=deltah+h_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltah=deltah+h_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltah=deltah+h_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltah=deltah+h_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltah=deltah+h_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltah=deltah+h_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltah=deltah+h_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltah=deltah+h_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltah=deltah+h_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltah=deltah+h_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltah=deltah+h_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltah=deltah+h_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltah=deltah+h_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltah=deltah+h_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltah=deltah+h_tt;
    }

    sug_enthalpy=deltah*1000;

    return sug_enthalpy;

}



double sug_entropy(int sequence_length, const std::string& specie)
{
    double sug_entropy;
    int i;
    double s_aa, s_ac, s_ag, s_at, s_ca, s_cc, s_cg, s_ct, s_ga, s_gc, s_gg, s_gt, s_ta, s_tc, s_tg, s_tt;
    double deltas;

    deltas=0;
    i=0;

    //Entropy
    //cal/(K mol)
    s_aa=s_tt=-21.9;
    s_at=-15.2;
    s_ta=-18.4;
    s_ca=s_tg=-21.0;
    s_gt=s_ac=-25.5;
    s_ct=s_ag=-16.4;
    s_ga=s_tc=-23.5;
    s_cg=-29.0;
    s_gc=-26.4;
    s_gg=s_cc=-28.4;


    //ENTROPY
    for(i=0; i<sequence_length-1; i++)
    {
        if (specie[i]=='A' && specie[i+1]=='A') deltas=deltas+s_aa;
        if (specie[i]=='A' && specie[i+1]=='C') deltas=deltas+s_ac;
        if (specie[i]=='A' && specie[i+1]=='G') deltas=deltas+s_ag;
        if (specie[i]=='A' && specie[i+1]=='T') deltas=deltas+s_at;

        if (specie[i]=='C' && specie[i+1]=='A') deltas=deltas+s_ca;
        if (specie[i]=='C' && specie[i+1]=='C') deltas=deltas+s_cc;
        if (specie[i]=='C' && specie[i+1]=='G') deltas=deltas+s_cg;
        if (specie[i]=='C' && specie[i+1]=='T') deltas=deltas+s_ct;

        if (specie[i]=='G' && specie[i+1]=='A') deltas=deltas+s_ga;
        if (specie[i]=='G' && specie[i+1]=='C') deltas=deltas+s_gc;
        if (specie[i]=='G' && specie[i+1]=='G') deltas=deltas+s_gg;
        if (specie[i]=='G' && specie[i+1]=='T') deltas=deltas+s_gt;

        if (specie[i]=='T' && specie[i+1]=='A') deltas=deltas+s_ta;
        if (specie[i]=='T' && specie[i+1]=='C') deltas=deltas+s_tc;
        if (specie[i]=='T' && specie[i+1]=='G') deltas=deltas+s_tg;
        if (specie[i]=='T' && specie[i+1]=='T') deltas=deltas+s_tt;

    }

    sug_entropy=deltas;

    return sug_entropy;

}



double sug_melting_curve(int sequence_length, const std::string& specie, double dna_conc)
{
    double f, keq, ctkeq;

    double deltah = sug_enthalpy(sequence_length, specie);
    double deltas = sug_entropy(sequence_length, specie);

    double t;
    double R=1.987; //cal/(K mol)

    t=0.0;
    std::ofstream fileout("sug_melting_curve.out");


    for (t=0.0; t<=700.0; t=t+0.5)
    {
        ctkeq=dna_conc*exp((deltas/R)-(deltah/(R*t)));
        f=(1+ctkeq-sqrt(1+2*ctkeq))/ctkeq;
        fileout << t << " " << f << std::endl;
    }
  
    return 0;

}

}

}
