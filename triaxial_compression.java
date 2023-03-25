/*
 * Copyright (c) 2022 Ehsan Haghighat & David Santillan
 * All rights reserved.
 *
 * triaxial_compression.java
 *
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

public class triaxial_compression {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.label("triaxial_compression.mph");

    model.param().set("G0b", "10 [MPa]");
    model.param().set("E0b", "2*G0b*(1+nu)");
    model.param().set("nu", "0.3");
    model.param().set("coh", "40 [kPa]");
    model.param().set("coh_r", "0 [kPa]");
    model.param().set("phi", "20 [deg]");
    model.param().set("phi_r", "20 [deg]");
    model.param().set("GII", "30 [J/m^2]");
    model.param().set("B", "80 [mm]");
    model.param().set("H", "170 [mm]");
    model.param().group().create("par2");
    model.param("par2").set("P0", "200 [kPa]");
    model.param("par2").set("SX0", "-P0");
    model.param("par2").set("SY0", "-P0*1.01");
    model.param("par2").set("SZ0", "nu*(SX0+SY0)");
    model.param().group().create("par3");
    model.param("par3").set("n", "2");
    model.param("par3").set("p", "1");
    model.param("par3").set("c0", "8/3");
    model.param("par3").set("K", "2*ell^2");
    model.param("par3").set("A", "1 [mm]");
    model.param("par3").set("ell", "2*A");
    model.param("par3").set("epsilon", "1e-6");
    model.param("par2").label("Initial stresses");
    model.param("par3").label("Phasefield ");

    model.component().create("comp1", false);

    model.component("comp1").geom().create("geom1", 2);

    model.component("comp1").curvedInterior(false);

    model.result().table().create("tbl1", "Table");
    model.result().table().create("evl2", "Table");

    model.component("comp1").func().create("int1", "Interpolation");
    model.component("comp1").func().create("an3", "Analytic");
    model.component("comp1").func().create("an4", "Analytic");
    model.component("comp1").func("int1").set("funcname", "Disp");
    model.component("comp1").func("int1").set("table", new String[][]{{"0", "0[mm]"}, {"10", "10[mm]"}, {"", ""}});
    model.component("comp1").func("int1").set("fununit", new String[]{"m"});
    model.component("comp1").func("int1").set("argunit", new String[]{"s"});
    model.component("comp1").func("an3").label("d_init");
    model.component("comp1").func("an3").set("funcname", "d_init");
    model.component("comp1").func("an3").set("expr", "exp(-2*((x-0.02)^2 + (y-0.04)^2)/0.001^2)");
    model.component("comp1").func("an3").set("args", new String[]{"x", "y"});
    model.component("comp1").func("an3").set("argunit", new String[]{"", ""});
    model.component("comp1").func("an3")
         .set("plotargs", new String[][]{{"x", "0", "80 [mm]"}, {"y", "0", "170 [mm]"}});
    model.component("comp1").func("an4").label("Analytic E0");
    model.component("comp1").func("an4").set("expr", "E0b*(1 - 0.2*exp(-((x-x0)^2 + (y-y0)^2)/0.000001))");
    model.component("comp1").func("an4").set("args", new String[]{"x", "y"});
    model.component("comp1").func("an4").set("argunit", new String[]{"", ""});
    model.component("comp1").func("an4")
         .set("plotargs", new String[][]{{"x", "0", "80 [mm]"}, {"y", "0", "170 [mm]"}});

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").repairTolType("relative");
    model.component("comp1").geom("geom1").create("r1", "Rectangle");
    model.component("comp1").geom("geom1").feature("r1").set("size", new String[]{"B", "H"});
    model.component("comp1").geom("geom1").create("pt1", "Point");
    model.component("comp1").geom("geom1").feature("pt1").set("p", new String[]{"B/2", "0"});
    model.component("comp1").geom("geom1").create("pt2", "Point");
    model.component("comp1").geom("geom1").feature("pt2").set("p", new String[]{"B/2", "H"});
    model.component("comp1").geom("geom1").feature("fin").set("repairtoltype", "relative");
    model.component("comp1").geom("geom1").run();

    model.component("comp1").variable().create("var6");
    model.component("comp1").variable("var6").set("x0", "B/2");
    model.component("comp1").variable("var6").set("y0", "H/2");
    model.component("comp1").variable("var6").set("E0", "E0b*(1 + 0.2*exp(-((x-x0)^2 + (y-y0)^2)/(1 [mm])^2))");
    model.component("comp1").variable("var6").set("K0", "E0/(3*(1-2*nu))");
    model.component("comp1").variable("var6").set("G0", "E0/(2*(1+nu))");
    model.component("comp1").variable("var6").set("L0", "K0 - 2*G0/3");
    model.component("comp1").variable("var6").set("U_imp", "Disp(t)");
    model.component("comp1").variable("var6").set("M", "GII/c0/ell");
    model.component("comp1").variable("var6").set("m", "M/Ht");
    model.component("comp1").variable().create("var5");
    model.component("comp1").variable("var5").set("e_vol", "solid.evol");
    model.component("comp1").variable("var5").set("ex", "solid.eXX - e_vol/3");
    model.component("comp1").variable("var5").set("ey", "solid.eYY - e_vol/3");
    model.component("comp1").variable("var5").set("ez", "solid.eZZ - e_vol/3");
    model.component("comp1").variable("var5").set("exy", "solid.eXY");
    model.component("comp1").variable("var5").set("e_eq", "sqrt(2/3*(ex^2 + ey^2 + ez^2 + 2*exy^2 + eps))");
    model.component("comp1").variable().create("var3");
    model.component("comp1").variable("var3").set("Sx", "K0*e_vol + 2*G0*ex + SX0");
    model.component("comp1").variable("var3").set("Sy", "K0*e_vol + 2*G0*ey + SY0");
    model.component("comp1").variable("var3").set("Sz", "K0*e_vol + 2*G0*ez + SZ0");
    model.component("comp1").variable("var3").set("sxy", "2*G0*exy");
    model.component("comp1").variable("var3").set("pm", "-(Sx+Sy+Sz)/3");
    model.component("comp1").variable("var3").set("sx", "Sx + pm");
    model.component("comp1").variable("var3").set("sy", "Sy + pm");
    model.component("comp1").variable("var3").set("sz", "Sz + pm");
    model.component("comp1").variable("var3").set("mises", "sqrt(1.5*(sx^2 + sy^2 + sz^2 + 2*sxy^2 + eps))");
    model.component("comp1").variable("var3").set("rinv", "4.5*(sx^3 + sy^3 + sz^3 + 3*sx*sxy^2 + 3*sxy^2*sy)");
    model.component("comp1").variable("var3").set("cos3th_trial", "rinv/(mises)^3");
    model.component("comp1").variable("var3")
         .set("cos3th", "if(abs(cos3th_trial)>0.999999999, sign(cos3th_trial), cos3th_trial)");
    model.component("comp1").variable("var3").set("Lth", "Lth_old", "acos(cos3th)/3");
    model.component("comp1").variable("var3")
         .set("Rmc", "sin(Lth + pi/3)/(sqrt(3)*cos(phi)) + cos(Lth + pi/3) * tan(phi) /3");
    model.component("comp1").variable("var3").set("ax", "sx/mises");
    model.component("comp1").variable("var3").set("ay", "sy/mises");
    model.component("comp1").variable("var3").set("axy", "sxy/mises");
    model.component("comp1").variable().create("var4");
    model.component("comp1").variable("var4").set("sigma_n", "pm");
    model.component("comp1").variable("var4").set("tau_b", "mises");
    model.component("comp1").variable("var4").set("tau_p", "(sigma_n*tan(phi) + coh)/Rmc");
    model.component("comp1").variable("var4").set("tau_r", "(sigma_n*tan(phi_r) + coh_r)/Rmc");
    model.component("comp1").variable("var4").set("f_mc", "tau_b - tau_p");
    model.component("comp1").variable("var4").set("gamma_t", "e_eq");
    model.component("comp1").variable("var4").set("Ht", "(tau_p - tau_r)^2/(6*G0)");
    model.component("comp1").variable("var4").set("Gamma", "Ht - (tau_p - tau_r)^2/(6*G0)");
    model.component("comp1").variable("var4").set("Hslip", "((tau_b-tau_r)^2 - (tau_p-tau_r)^2)/(6*G0)");
    model.component("comp1").variable("var4").set("Henergy", "Hslip");
    model.component("comp1").variable("var4").set("H_max", "max(0, max(H_max_old, Henergy))");
    model.component("comp1").variable("var4").set("gd", "(1-d)^n/((1-d)^n + m*d*(1+p*d))");
    model.component("comp1").variable("var4").set("dgd", "d(gd, d)");
    model.component("comp1").variable("var4").set("wd", "d");
    model.component("comp1").variable("var4").set("dwd", "1");

    model.view().create("view2", 3);
    model.view().create("view3", 2);
    model.view().create("view4", 3);

    model.component("comp1").physics().create("solid", "SolidMechanics", "geom1");
    model.component("comp1").physics("solid").feature("lemm1").create("iss1", "InitialStressandStrain", 2);
    model.component("comp1").physics("solid").feature("lemm1").create("exs1", "ExternalStress", 2);
    model.component("comp1").physics("solid").create("disp8", "Displacement1", 1);
    model.component("comp1").physics("solid").feature("disp8").selection().set(2, 4);
    model.component("comp1").physics("solid").create("disp7", "Displacement1", 1);
    model.component("comp1").physics("solid").feature("disp7").selection().set(3, 5);
    model.component("comp1").physics("solid").create("disp10", "Displacement0", 0);
    model.component("comp1").physics("solid").feature("disp10").selection().set(4);
    model.component("comp1").physics("solid").create("bndl1", "BoundaryLoad", 1);
    model.component("comp1").physics("solid").feature("bndl1").selection().set(1, 6);
    model.component("comp1").physics().create("c", "CoefficientFormPDE", "geom1");
    model.component("comp1").physics("c").field("dimensionless").field("d");
    model.component("comp1").physics("c").field("dimensionless").component(new String[]{"d"});
    model.component("comp1").physics("c").create("dir1", "DirichletBoundary", 1);
    model.component("comp1").physics().create("dode", "DomainODE", "geom1");
    model.component("comp1").physics("dode").field("dimensionless").field("H_max_old");
    model.component("comp1").physics("dode").field("dimensionless").component(new String[]{"H_max_old"});
    model.component("comp1").physics().create("dode2", "DomainODE", "geom1");
    model.component("comp1").physics("dode2").field("dimensionless").field("Lth_old");
    model.component("comp1").physics("dode2").field("dimensionless").component(new String[]{"Lth_old"});

    model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh1").feature("ftri1").selection().geom("geom1");

    model.component("comp1").probe().create("bnd1", "Boundary");
    model.component("comp1").probe().create("bnd3", "Boundary");
    model.component("comp1").probe().create("bnd2", "Boundary");
    model.component("comp1").probe().create("bnd4", "Boundary");
    model.component("comp1").probe().create("dom1", "Domain");
    model.component("comp1").probe().create("dom3", "Domain");
    model.component("comp1").probe().create("dom4", "Domain");
    model.component("comp1").probe().create("dom5", "Domain");
    model.component("comp1").probe().create("dom2", "Domain");
    model.component("comp1").probe("bnd1").selection().set(2, 4);
    model.component("comp1").probe("bnd3").selection().set(2, 4);
    model.component("comp1").probe("bnd2").selection().set(3, 5);
    model.component("comp1").probe("bnd4").selection().set(3, 5);

    model.result().table("tbl1").label("Probe Table 1");
    model.result().table("evl2").label("Evaluation 2D");
    model.result().table("evl2").comments("Interactive 2D values");

    model.thermodynamics().label("Thermodynamics Package");

    model.component("comp1").variable("var6").label("Params");
    model.component("comp1").variable("var5").label("Strain Components");
    model.component("comp1").variable("var3").label("Stress Components");
    model.component("comp1").variable("var4").label("Phase Field Global");

    model.component("comp1").view("view1").axis().set("xmin", -0.12573018670082092);
    model.component("comp1").view("view1").axis().set("xmax", 0.20573018491268158);
    model.component("comp1").view("view1").axis().set("ymin", -0.004249997437000275);
    model.component("comp1").view("view1").axis().set("ymax", 0.17425000667572021);
    model.view("view3").axis().set("xmin", -2.115384578704834);
    model.view("view3").axis().set("xmax", 2.115384578704834);

    model.component("comp1").physics("solid").prop("ShapeProperty").set("order_displacement", 1);
    model.component("comp1").physics("solid").prop("ShapeProperty").set("valueType", "real");
    model.component("comp1").physics("solid").prop("d").set("d", "B");
    model.component("comp1").physics("solid").prop("StructuralTransientBehavior")
         .set("StructuralTransientBehavior", "Quasistatic");
    model.component("comp1").physics("solid").feature("lemm1").set("D_mat", "userdef");
    model.component("comp1").physics("solid").feature("lemm1")
         .set("D", new String[][]{{"2*G0+L0"}, 
         {"L0"}, 
         {"L0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"L0"}, 
         {"2*G0+L0"}, 
         {"L0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"L0"}, 
         {"L0"}, 
         {"2*G0+L0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"2*G0*g(d)"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"2*G0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"0"}, 
         {"2*G0"}});
    model.component("comp1").physics("solid").feature("lemm1").set("IsotropicOption", "KG");
    model.component("comp1").physics("solid").feature("lemm1").set("E_mat", "userdef");
    model.component("comp1").physics("solid").feature("lemm1").set("E", "E0");
    model.component("comp1").physics("solid").feature("lemm1").set("nu_mat", "userdef");
    model.component("comp1").physics("solid").feature("lemm1").set("nu", "nu");
    model.component("comp1").physics("solid").feature("lemm1").set("K_mat", "userdef");
    model.component("comp1").physics("solid").feature("lemm1").set("K", "K0");
    model.component("comp1").physics("solid").feature("lemm1").set("G_mat", "userdef");
    model.component("comp1").physics("solid").feature("lemm1").set("G", "G0*gd");
    model.component("comp1").physics("solid").feature("lemm1").set("lambLame_mat", "userdef");
    model.component("comp1").physics("solid").feature("lemm1").set("lambLame", "L0");
    model.component("comp1").physics("solid").feature("lemm1").set("muLame_mat", "userdef");
    model.component("comp1").physics("solid").feature("lemm1").set("muLame", "G0*(1-d)^2");
    model.component("comp1").physics("solid").feature("lemm1").set("rho_mat", "userdef");
    model.component("comp1").physics("solid").feature("lemm1").set("minput_temperature_src", "userdef");
    model.component("comp1").physics("solid").feature("lemm1")
         .set("minput_strainreferencetemperature_src", "userdef");
    model.component("comp1").physics("solid").feature("lemm1").feature("iss1")
         .set("Sil", new String[][]{{"SX0"}, {"0"}, {"0"}, {"0"}, {"SY0"}, {"0"}, {"0"}, {"0"}, {"SZ0"}});
    model.component("comp1").physics("solid").feature("lemm1").feature("exs1")
         .set("Sext", new String[][]{{"(1-gd)*tau_r*ax"}, {"(1-gd)*tau_r*axy"}, {"0"}, {"(1-gd)*tau_r*axy"}, {"(1-gd)*tau_r*ay"}, {"0"}, {"0"}, {"0"}, {"0"}});
    model.component("comp1").physics("solid").feature("dcnt1").set("pairDisconnect", true);
    model.component("comp1").physics("solid").feature("dcnt1").set("integrationOrder", 2);
    model.component("comp1").physics("solid").feature("dcnt1").label("Contact");
    model.component("comp1").physics("solid").feature("dcont1").set("pairDisconnect", true);
    model.component("comp1").physics("solid").feature("dcont1").label("Continuity");
    model.component("comp1").physics("solid").feature("disp8").set("Direction", new int[][]{{0}, {1}, {0}});
    model.component("comp1").physics("solid").feature("disp8").label("Prescribed Displacement Bot");
    model.component("comp1").physics("solid").feature("disp7").set("Direction", new int[][]{{0}, {1}, {0}});
    model.component("comp1").physics("solid").feature("disp7").set("U0", new String[][]{{"0"}, {"-Disp(t)"}, {"0"}});
    model.component("comp1").physics("solid").feature("disp7").label("Prescribed Displacement Top");
    model.component("comp1").physics("solid").feature("disp10").set("Direction", new int[][]{{1}, {0}, {0}});
    model.component("comp1").physics("solid").feature("bndl1").set("LoadType", "FollowerPressure");
    model.component("comp1").physics("solid").feature("bndl1").set("FollowerPressure", "P0");
    model.component("comp1").physics("solid").feature("bndl1").label("Boundary Confining");
    model.component("comp1").physics("c").prop("ShapeProperty").set("order", 1);
    model.component("comp1").physics("c").prop("ShapeProperty").set("valueType", "real");
    model.component("comp1").physics("c").feature("cfeq1").set("c", new String[][]{{"M*K", "0", "0", "M*K"}});
    model.component("comp1").physics("c").feature("cfeq1").set("f", "-M*dwd - dgd*(Ht + H_max)");
    model.component("comp1").physics("c").feature("dcont1").set("pairDisconnect", true);
    model.component("comp1").physics("c").feature("dcont1").label("Continuity");
    model.component("comp1").physics("c").feature("dir1").set("r", 1);
    model.component("comp1").physics("c").feature("dir1").active(false);
    model.component("comp1").physics("dode").prop("ShapeProperty").set("shapeFunctionType", "shgp");
    model.component("comp1").physics("dode").prop("ShapeProperty").set("order", 0);
    model.component("comp1").physics("dode").prop("Units").set("SourceTermQuantity", "dimensionless");
    model.component("comp1").physics("dode").feature("dode1")
         .set("f", "H_max_old-nojac(if(Henergy>H_max_old, Henergy,H_max_old))");
    model.component("comp1").physics("dode").feature("dode1").set("da", 0);
    model.component("comp1").physics("dode2").prop("ShapeProperty").set("shapeFunctionType", "shgp");
    model.component("comp1").physics("dode2").prop("ShapeProperty").set("order", 0);
    model.component("comp1").physics("dode2").prop("Units").set("CustomSourceTermUnit", 1);
    model.component("comp1").physics("dode2").feature("dode1").set("f", "Lth_old - nojac(acos(cos3th)/3)");
    model.component("comp1").physics("dode2").feature("dode1").set("da", 0);

    model.component("comp1").mesh("mesh1").feature("size").set("hauto", 1);
    model.component("comp1").mesh("mesh1").feature("ftri1").set("smoothmaxiter", 8);
    model.component("comp1").mesh("mesh1").feature("ftri1").set("smoothmaxdepth", 8);
    model.component("comp1").mesh("mesh1").run();

    model.component("comp1").probe("bnd1").label("Boundary Probe RFy");
    model.component("comp1").probe("bnd1").set("type", "integral");
    model.component("comp1").probe("bnd1").set("probename", "bnd_rfy");
    model.component("comp1").probe("bnd1").set("expr", "solid.RFy/B");
    model.component("comp1").probe("bnd1").set("unit", "N/m");
    model.component("comp1").probe("bnd1").set("descr", "solid.RFy/B");
    model.component("comp1").probe("bnd1").set("method", "summation");
    model.component("comp1").probe("bnd1").set("frame", "geometry");
    model.component("comp1").probe("bnd1").set("table", "tbl1");
    model.component("comp1").probe("bnd1").set("window", "window2");
    model.component("comp1").probe("bnd3").label("Boundary Probe RFx");
    model.component("comp1").probe("bnd3").set("type", "integral");
    model.component("comp1").probe("bnd3").set("probename", "bnd_rfx");
    model.component("comp1").probe("bnd3").set("expr", "solid.RFx");
    model.component("comp1").probe("bnd3").set("unit", "N");
    model.component("comp1").probe("bnd3").set("descr", "Reaction force, x component");
    model.component("comp1").probe("bnd3").set("method", "summation");
    model.component("comp1").probe("bnd3").set("frame", "geometry");
    model.component("comp1").probe("bnd3").set("table", "tbl1");
    model.component("comp1").probe("bnd3").set("window", "window2");
    model.component("comp1").probe("bnd2").label("Boundary Probe UX");
    model.component("comp1").probe("bnd2").set("probename", "bnd_ux");
    model.component("comp1").probe("bnd2").set("expr", "u");
    model.component("comp1").probe("bnd2").set("unit", "1");
    model.component("comp1").probe("bnd2").set("descr", "Displacement field, X component");
    model.component("comp1").probe("bnd2").set("method", "summation");
    model.component("comp1").probe("bnd2").set("frame", "geometry");
    model.component("comp1").probe("bnd2").set("table", "tbl1");
    model.component("comp1").probe("bnd2").set("window", "window2");
    model.component("comp1").probe("bnd4").label("Boundary Probe UY");
    model.component("comp1").probe("bnd4").set("probename", "bnd_uy");
    model.component("comp1").probe("bnd4").set("expr", "v");
    model.component("comp1").probe("bnd4").set("unit", "1");
    model.component("comp1").probe("bnd4").set("descr", "Displacement field, Y component");
    model.component("comp1").probe("bnd4").set("method", "summation");
    model.component("comp1").probe("bnd4").set("frame", "geometry");
    model.component("comp1").probe("bnd4").set("table", "tbl1");
    model.component("comp1").probe("bnd4").set("window", "window2");
    model.component("comp1").probe("dom1").label("Domain Probe Sy");
    model.component("comp1").probe("dom1").set("probename", "sy_mean");
    model.component("comp1").probe("dom1").set("expr", "solid.sy");
    model.component("comp1").probe("dom1").set("unit", "N/m^2");
    model.component("comp1").probe("dom1").set("descr", "Stress tensor, y component");
    model.component("comp1").probe("dom1").set("table", "tbl1");
    model.component("comp1").probe("dom1").set("window", "window2");
    model.component("comp1").probe("dom3").label("Domain Probe Sx");
    model.component("comp1").probe("dom3").set("probename", "sx_mean");
    model.component("comp1").probe("dom3").set("expr", "solid.sx");
    model.component("comp1").probe("dom3").set("unit", "N/m^2");
    model.component("comp1").probe("dom3").set("descr", "Stress tensor, x component");
    model.component("comp1").probe("dom3").set("table", "tbl1");
    model.component("comp1").probe("dom3").set("window", "window2");
    model.component("comp1").probe("dom4").label("Domain Probe Sxy");
    model.component("comp1").probe("dom4").set("probename", "sxy_mean");
    model.component("comp1").probe("dom4").set("expr", "solid.sxy");
    model.component("comp1").probe("dom4").set("unit", "N/m^2");
    model.component("comp1").probe("dom4").set("descr", "Stress tensor, xy component");
    model.component("comp1").probe("dom4").set("table", "tbl1");
    model.component("comp1").probe("dom4").set("window", "window2");
    model.component("comp1").probe("dom5").label("Domain Probe Mises");
    model.component("comp1").probe("dom5").set("probename", "mises_mean");
    model.component("comp1").probe("dom5").set("expr", "solid.mises");
    model.component("comp1").probe("dom5").set("unit", "N/m^2");
    model.component("comp1").probe("dom5").set("descr", "von Mises stress");
    model.component("comp1").probe("dom5").set("table", "tbl1");
    model.component("comp1").probe("dom5").set("window", "window2");
    model.component("comp1").probe("dom2").label("Domain Probe Gamma");
    model.component("comp1").probe("dom2").set("probename", "gamma_mean");
    model.component("comp1").probe("dom2").set("expr", "solid.eXY");
    model.component("comp1").probe("dom2").set("unit", "1");
    model.component("comp1").probe("dom2").set("descr", "Strain tensor, XY component");
    model.component("comp1").probe("dom2").set("table", "tbl1");
    model.component("comp1").probe("dom2").set("window", "window2");

    model.study().create("std2");
    model.study("std2").create("time", "Transient");

    model.sol().create("sol1");
    model.sol("sol1").study("std2");
    model.sol("sol1").attach("std2");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").create("t1", "Time");
    model.sol("sol1").feature("t1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("t1").create("se1", "Segregated");
    model.sol("sol1").feature("t1").create("ps1", "PreviousSolution");
    model.sol("sol1").feature("t1").feature("se1").create("c2", "SegregatedStep");
    model.sol("sol1").feature("t1").feature().remove("fcDef");

    model.result().dataset().remove("dset1");
    model.result().create("pg5", "PlotGroup1D");
    model.result("pg5").set("probetag", "window2");
    model.result("pg5").create("tblp1", "Table");
    model.result("pg5").feature("tblp1").set("probetag", "bnd1,bnd3,bnd2,bnd4,dom1,dom3,dom4,dom2,dom5");

    model.component("comp1").probe("bnd1").genResult(null);
    model.component("comp1").probe("bnd3").genResult(null);
    model.component("comp1").probe("bnd2").genResult(null);
    model.component("comp1").probe("bnd4").genResult(null);
    model.component("comp1").probe("dom1").genResult(null);
    model.component("comp1").probe("dom3").genResult(null);
    model.component("comp1").probe("dom4").genResult(null);
    model.component("comp1").probe("dom5").genResult(null);
    model.component("comp1").probe("dom2").genResult(null);

    model.study("std2").label("Study 1");
    model.study("std2").feature("time").set("tlist", "range(0,.1,10)");
    model.study("std2").feature("time").set("plot", true);
    model.study("std2").feature("time").set("plotfreq", "tsteps");
    model.study("std2").feature("time").set("useinitsol", true);
    model.study("std2").feature("time").set("initmethod", "sol");
    model.study("std2").feature("time").set("usesol", true);
    model.study("std2").feature("time").set("notsolmethod", "sol");

    model.sol("sol1").attach("std2");
    model.sol("sol1").feature("st1").label("Compile Equations: Time Dependent");
    model.sol("sol1").feature("v1").label("Dependent Variables 1.1");
    model.sol("sol1").feature("v1").set("initmethod", "sol");
    model.sol("sol1").feature("v1").set("resscalemethod", "auto");
    model.sol("sol1").feature("v1").set("notsolmethod", "sol");
    model.sol("sol1").feature("v1").set("clist", new String[]{"range(0,.1,10)", "0.01[s]"});
    model.sol("sol1").feature("v1").feature("comp1_u").set("scalemethod", "manual");
    model.sol("sol1").feature("v1").feature("comp1_u").set("scaleval", "1e-2*1.8027756377319948");
    model.sol("sol1").feature("t1").label("Time-Dependent Solver 1.1");
    model.sol("sol1").feature("t1").set("control", "user");
    model.sol("sol1").feature("t1").set("tlist", "range(0,.1,10)");
    model.sol("sol1").feature("t1").set("tout", "tsteps");
    model.sol("sol1").feature("t1").set("atolglobalvaluemethod", "manual");
    model.sol("sol1").feature("t1")
         .set("atolvaluemethod", new String[]{"comp1_u", "manual", "comp1_H_max_old", "manual", "comp1_Lth_old", "factor", "comp1_d", "factor"});
    model.sol("sol1").feature("t1").set("maxstepconstraintbdf", "const");
    model.sol("sol1").feature("t1").set("maxstepbdf", 2);
    model.sol("sol1").feature("t1").set("plot", true);
    model.sol("sol1").feature("t1").set("plotfreq", "tsteps");
    model.sol("sol1").feature("t1").feature("dDef").label("Direct 1");
    model.sol("sol1").feature("t1").feature("dDef").set("linsolver", "pardiso");
    model.sol("sol1").feature("t1").feature("dDef").set("ooc", false);
    model.sol("sol1").feature("t1").feature("dDef").set("pardreorder", "ndmt");
    model.sol("sol1").feature("t1").feature("dDef").set("pardmtsolve", false);
    model.sol("sol1").feature("t1").feature("dDef").set("rhob", 400);
    model.sol("sol1").feature("t1").feature("aDef").label("Advanced 1");
    model.sol("sol1").feature("t1").feature("fc1").label("Fully Coupled 1.1");
    model.sol("sol1").feature("t1").feature("se1").label("Segregated 1.1");
    model.sol("sol1").feature("t1").feature("se1").set("segterm", "itertol");
    model.sol("sol1").feature("t1").feature("se1").set("segiter", 15);
    model.sol("sol1").feature("t1").feature("se1").set("plot", true);
    model.sol("sol1").feature("t1").feature("se1").feature("ssDef").label("Segregated Step 1");
    model.sol("sol1").feature("t1").feature("se1").feature("ssDef").set("segvar", new String[]{"comp1_u"});
    model.sol("sol1").feature("t1").feature("se1").feature("ssDef").set("subdtech", "auto");
    model.sol("sol1").feature("t1").feature("se1").feature("ssDef").set("maxsubiter", 50);
    model.sol("sol1").feature("t1").feature("se1").feature("c2").label("Segregated Step 2.1");
    model.sol("sol1").feature("t1").feature("se1").feature("c2").set("segvar", new String[]{"comp1_d"});
    model.sol("sol1").feature("t1").feature("ps1").label("Previous Solution 1.1");
    model.sol("sol1").feature("t1").feature("ps1").set("prevcomp", new String[]{"comp1_H_max_old", "comp1_Lth_old"});
    model.sol("sol1").runAll();

    model.result().dataset().remove("dset1");
    model.result().dataset().remove("dset2");
    model.result("pg5").label("Probe Plot Group 5");
    model.result("pg5").set("solrepresentation", "solnum");
    model.result("pg5").set("xlabel", "Time (s)");
    model.result("pg5").set("ylabel", "Dependent variable d (N/m), Boundary Probe RFy");
    model.result("pg5").set("axislimits", true);
    model.result("pg5").set("xmin", 0);
    model.result("pg5").set("xmax", 10);
    model.result("pg5").set("ymin", 0);
    model.result("pg5").set("ymax", 23000);
    model.result("pg5").set("windowtitle", "Probe Plot 2");
    model.result("pg5").set("xlabelactive", false);
    model.result("pg5").set("ylabelactive", false);
    model.result("pg5").feature("tblp1").label("Probe Table Graph 1");
    model.result("pg5").feature("tblp1").set("legend", true);

    return model;
  }

  public static void main(String[] args) {
    run();
  }

}
