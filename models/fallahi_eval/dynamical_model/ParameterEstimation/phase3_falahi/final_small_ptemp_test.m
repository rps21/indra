function [err, timepoints, species_out, observables_out ] = final_small_ptemp_test( timepoints, species_init, parameters, suppress_plot )
%FINAL_SMALL_PTEMP_TEST Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'final_small_ptemp_test' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the CVode library interfaced
%   to MATLAB via the MEX interface. Before running this script, the model
%   source in file final_small_ptemp_test_cvode.c must be compiled (see that file for details).
%   FINAL_SMALL_PTEMP_TEST returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = final_small_ptemp_test( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   timepoints      : column vector of time points returned by integrator.
%   species_init    : row vector of 725 initial species populations.
%   parameters      : row vector of 77 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 2 )
    species_init = [];
end

if ( nargin < 3 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 0;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 1.000000e+02, 1e-6, 1e-1, 1e2, 1.000000e-06, 1.000000e-01, 1.000000e-06, 1.000000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.000000e+04, 1.0e6, 1, 1, 1, 1, 1 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 77  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 77].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 725  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 725].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,10,20+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'kcat_braf_act', 'kf_mtorc2_akt', 'kr_mtorc2_akt', 'kphos_s473', 'kf_egfr_grb2_bind', 'kr_egfr_grb2_unbind', 'kf_egf_egfr_bind', 'kr_egf_egfr_bind', 'kf_egfr_akt', 'kr_egfr_akt', 'kc_egfr_akt', 'kf_akt_raf1', 'kr_akt_raf1', 'kc_akt_raf1', 'kf_raf1_map2k1', 'kr_raf1_map2k1', 'kc__raf1_map2k1', 'kt_raf1', 'kf_map2k1_mapk3', 'kc_map2k1_mapk3', 'kr_map2k1_mapk3', 'kf_mapk1_map2k1', 'kr_mapk1_map2k1', 'kc_mapk1_map2k1', 'kf_mapk1_mapk3', 'kr_mapk1_mapk3', 'kf_mtor_rictor', 'kr_mtor_rictor', 'kf_mtor_rptor', 'kr_mtor_rptor', 'kf_egfr_rictor', 'kc_egfr_rictor', 'kr__egfr_rictor', 'kf_mtor_s6k', 'kc__mtor_s6k', 'kr__mtor_s6k', 'kf_rictor_s6k', 'kr_rictor_s6k', 'kc_rictor_s6k', 'kf_s6k_s6', 'kc_s6k_s6', 'kr_s6k_s6', 'kf_akt_cdkn1b', 'kc_akt_cdkn1b', 'kr_akt_cdkn1b', 'kf_grb2_sos', 'kr_grb2_sos', 'kf_kras_sos', 'kr__kras_sos', 'kf_kras_braf', 'kr_kras_braf', 'kc_kras_braf', 'kf_braf_map2k1', 'kr_braf_map2k1', 'kc_braf_map2k1', 'RPTOR_0', 'RAF1_0', 'GRB2_0', 'MAPK1_0', 'SOS1_0', 'MTOR_0', 'RPS6KB1_0', 'MAPK3_0', 'BRAF_0', 'AKT1_0', 'EGFR_0', 'KRAS_0', 'RICTOR_0', 'MAP2K1_0', 'CDKN1B_0', 'RPS6_0', 'EGF_0', 'drug_0', 'kf_drug_braf', 'kr_drug_braf', 'kf_drug_raf1', 'kr_drug_raf1' };



%% Integrate Network Model
try 
    % run simulation
    [err, species_out, observables_out] = final_small_ptemp_test_cvode( timepoints, species_init, parameters );
catch
    fprintf( 1, 'Error: some problem integrating ODE network! (CVODE exitflag %d)\n', err );
    err = 1;
    return;
end



%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'pmek', 'perk', 'pakt308', 'pakt473', 'pmtors2448', 'pp70S6K', 'ps6', 'p27' };

    % construct figure
    plot(timepoints,observables_out);
    title('final_small_ptemp_test observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end



%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%



% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,725);
    species_init(1) = params(72);
    species_init(2) = params(56);
    species_init(3) = params(57);
    species_init(4) = params(58);
    species_init(5) = params(59);
    species_init(6) = params(60);
    species_init(7) = params(61);
    species_init(8) = params(62);
    species_init(9) = params(63);
    species_init(10) = params(64);
    species_init(11) = params(65);
    species_init(12) = params(66);
    species_init(13) = params(67);
    species_init(14) = params(68);
    species_init(15) = params(69);
    species_init(16) = params(70);
    species_init(17) = params(71);
    species_init(18) = params(73);
    species_init(19) = 0;
    species_init(20) = 0;
    species_init(21) = 0;
    species_init(22) = 0;
    species_init(23) = 0;
    species_init(24) = 0;
    species_init(25) = 0;
    species_init(26) = 0;
    species_init(27) = 0;
    species_init(28) = 0;
    species_init(29) = 0;
    species_init(30) = 0;
    species_init(31) = 0;
    species_init(32) = 0;
    species_init(33) = 0;
    species_init(34) = 0;
    species_init(35) = 0;
    species_init(36) = 0;
    species_init(37) = 0;
    species_init(38) = 0;
    species_init(39) = 0;
    species_init(40) = 0;
    species_init(41) = 0;
    species_init(42) = 0;
    species_init(43) = 0;
    species_init(44) = 0;
    species_init(45) = 0;
    species_init(46) = 0;
    species_init(47) = 0;
    species_init(48) = 0;
    species_init(49) = 0;
    species_init(50) = 0;
    species_init(51) = 0;
    species_init(52) = 0;
    species_init(53) = 0;
    species_init(54) = 0;
    species_init(55) = 0;
    species_init(56) = 0;
    species_init(57) = 0;
    species_init(58) = 0;
    species_init(59) = 0;
    species_init(60) = 0;
    species_init(61) = 0;
    species_init(62) = 0;
    species_init(63) = 0;
    species_init(64) = 0;
    species_init(65) = 0;
    species_init(66) = 0;
    species_init(67) = 0;
    species_init(68) = 0;
    species_init(69) = 0;
    species_init(70) = 0;
    species_init(71) = 0;
    species_init(72) = 0;
    species_init(73) = 0;
    species_init(74) = 0;
    species_init(75) = 0;
    species_init(76) = 0;
    species_init(77) = 0;
    species_init(78) = 0;
    species_init(79) = 0;
    species_init(80) = 0;
    species_init(81) = 0;
    species_init(82) = 0;
    species_init(83) = 0;
    species_init(84) = 0;
    species_init(85) = 0;
    species_init(86) = 0;
    species_init(87) = 0;
    species_init(88) = 0;
    species_init(89) = 0;
    species_init(90) = 0;
    species_init(91) = 0;
    species_init(92) = 0;
    species_init(93) = 0;
    species_init(94) = 0;
    species_init(95) = 0;
    species_init(96) = 0;
    species_init(97) = 0;
    species_init(98) = 0;
    species_init(99) = 0;
    species_init(100) = 0;
    species_init(101) = 0;
    species_init(102) = 0;
    species_init(103) = 0;
    species_init(104) = 0;
    species_init(105) = 0;
    species_init(106) = 0;
    species_init(107) = 0;
    species_init(108) = 0;
    species_init(109) = 0;
    species_init(110) = 0;
    species_init(111) = 0;
    species_init(112) = 0;
    species_init(113) = 0;
    species_init(114) = 0;
    species_init(115) = 0;
    species_init(116) = 0;
    species_init(117) = 0;
    species_init(118) = 0;
    species_init(119) = 0;
    species_init(120) = 0;
    species_init(121) = 0;
    species_init(122) = 0;
    species_init(123) = 0;
    species_init(124) = 0;
    species_init(125) = 0;
    species_init(126) = 0;
    species_init(127) = 0;
    species_init(128) = 0;
    species_init(129) = 0;
    species_init(130) = 0;
    species_init(131) = 0;
    species_init(132) = 0;
    species_init(133) = 0;
    species_init(134) = 0;
    species_init(135) = 0;
    species_init(136) = 0;
    species_init(137) = 0;
    species_init(138) = 0;
    species_init(139) = 0;
    species_init(140) = 0;
    species_init(141) = 0;
    species_init(142) = 0;
    species_init(143) = 0;
    species_init(144) = 0;
    species_init(145) = 0;
    species_init(146) = 0;
    species_init(147) = 0;
    species_init(148) = 0;
    species_init(149) = 0;
    species_init(150) = 0;
    species_init(151) = 0;
    species_init(152) = 0;
    species_init(153) = 0;
    species_init(154) = 0;
    species_init(155) = 0;
    species_init(156) = 0;
    species_init(157) = 0;
    species_init(158) = 0;
    species_init(159) = 0;
    species_init(160) = 0;
    species_init(161) = 0;
    species_init(162) = 0;
    species_init(163) = 0;
    species_init(164) = 0;
    species_init(165) = 0;
    species_init(166) = 0;
    species_init(167) = 0;
    species_init(168) = 0;
    species_init(169) = 0;
    species_init(170) = 0;
    species_init(171) = 0;
    species_init(172) = 0;
    species_init(173) = 0;
    species_init(174) = 0;
    species_init(175) = 0;
    species_init(176) = 0;
    species_init(177) = 0;
    species_init(178) = 0;
    species_init(179) = 0;
    species_init(180) = 0;
    species_init(181) = 0;
    species_init(182) = 0;
    species_init(183) = 0;
    species_init(184) = 0;
    species_init(185) = 0;
    species_init(186) = 0;
    species_init(187) = 0;
    species_init(188) = 0;
    species_init(189) = 0;
    species_init(190) = 0;
    species_init(191) = 0;
    species_init(192) = 0;
    species_init(193) = 0;
    species_init(194) = 0;
    species_init(195) = 0;
    species_init(196) = 0;
    species_init(197) = 0;
    species_init(198) = 0;
    species_init(199) = 0;
    species_init(200) = 0;
    species_init(201) = 0;
    species_init(202) = 0;
    species_init(203) = 0;
    species_init(204) = 0;
    species_init(205) = 0;
    species_init(206) = 0;
    species_init(207) = 0;
    species_init(208) = 0;
    species_init(209) = 0;
    species_init(210) = 0;
    species_init(211) = 0;
    species_init(212) = 0;
    species_init(213) = 0;
    species_init(214) = 0;
    species_init(215) = 0;
    species_init(216) = 0;
    species_init(217) = 0;
    species_init(218) = 0;
    species_init(219) = 0;
    species_init(220) = 0;
    species_init(221) = 0;
    species_init(222) = 0;
    species_init(223) = 0;
    species_init(224) = 0;
    species_init(225) = 0;
    species_init(226) = 0;
    species_init(227) = 0;
    species_init(228) = 0;
    species_init(229) = 0;
    species_init(230) = 0;
    species_init(231) = 0;
    species_init(232) = 0;
    species_init(233) = 0;
    species_init(234) = 0;
    species_init(235) = 0;
    species_init(236) = 0;
    species_init(237) = 0;
    species_init(238) = 0;
    species_init(239) = 0;
    species_init(240) = 0;
    species_init(241) = 0;
    species_init(242) = 0;
    species_init(243) = 0;
    species_init(244) = 0;
    species_init(245) = 0;
    species_init(246) = 0;
    species_init(247) = 0;
    species_init(248) = 0;
    species_init(249) = 0;
    species_init(250) = 0;
    species_init(251) = 0;
    species_init(252) = 0;
    species_init(253) = 0;
    species_init(254) = 0;
    species_init(255) = 0;
    species_init(256) = 0;
    species_init(257) = 0;
    species_init(258) = 0;
    species_init(259) = 0;
    species_init(260) = 0;
    species_init(261) = 0;
    species_init(262) = 0;
    species_init(263) = 0;
    species_init(264) = 0;
    species_init(265) = 0;
    species_init(266) = 0;
    species_init(267) = 0;
    species_init(268) = 0;
    species_init(269) = 0;
    species_init(270) = 0;
    species_init(271) = 0;
    species_init(272) = 0;
    species_init(273) = 0;
    species_init(274) = 0;
    species_init(275) = 0;
    species_init(276) = 0;
    species_init(277) = 0;
    species_init(278) = 0;
    species_init(279) = 0;
    species_init(280) = 0;
    species_init(281) = 0;
    species_init(282) = 0;
    species_init(283) = 0;
    species_init(284) = 0;
    species_init(285) = 0;
    species_init(286) = 0;
    species_init(287) = 0;
    species_init(288) = 0;
    species_init(289) = 0;
    species_init(290) = 0;
    species_init(291) = 0;
    species_init(292) = 0;
    species_init(293) = 0;
    species_init(294) = 0;
    species_init(295) = 0;
    species_init(296) = 0;
    species_init(297) = 0;
    species_init(298) = 0;
    species_init(299) = 0;
    species_init(300) = 0;
    species_init(301) = 0;
    species_init(302) = 0;
    species_init(303) = 0;
    species_init(304) = 0;
    species_init(305) = 0;
    species_init(306) = 0;
    species_init(307) = 0;
    species_init(308) = 0;
    species_init(309) = 0;
    species_init(310) = 0;
    species_init(311) = 0;
    species_init(312) = 0;
    species_init(313) = 0;
    species_init(314) = 0;
    species_init(315) = 0;
    species_init(316) = 0;
    species_init(317) = 0;
    species_init(318) = 0;
    species_init(319) = 0;
    species_init(320) = 0;
    species_init(321) = 0;
    species_init(322) = 0;
    species_init(323) = 0;
    species_init(324) = 0;
    species_init(325) = 0;
    species_init(326) = 0;
    species_init(327) = 0;
    species_init(328) = 0;
    species_init(329) = 0;
    species_init(330) = 0;
    species_init(331) = 0;
    species_init(332) = 0;
    species_init(333) = 0;
    species_init(334) = 0;
    species_init(335) = 0;
    species_init(336) = 0;
    species_init(337) = 0;
    species_init(338) = 0;
    species_init(339) = 0;
    species_init(340) = 0;
    species_init(341) = 0;
    species_init(342) = 0;
    species_init(343) = 0;
    species_init(344) = 0;
    species_init(345) = 0;
    species_init(346) = 0;
    species_init(347) = 0;
    species_init(348) = 0;
    species_init(349) = 0;
    species_init(350) = 0;
    species_init(351) = 0;
    species_init(352) = 0;
    species_init(353) = 0;
    species_init(354) = 0;
    species_init(355) = 0;
    species_init(356) = 0;
    species_init(357) = 0;
    species_init(358) = 0;
    species_init(359) = 0;
    species_init(360) = 0;
    species_init(361) = 0;
    species_init(362) = 0;
    species_init(363) = 0;
    species_init(364) = 0;
    species_init(365) = 0;
    species_init(366) = 0;
    species_init(367) = 0;
    species_init(368) = 0;
    species_init(369) = 0;
    species_init(370) = 0;
    species_init(371) = 0;
    species_init(372) = 0;
    species_init(373) = 0;
    species_init(374) = 0;
    species_init(375) = 0;
    species_init(376) = 0;
    species_init(377) = 0;
    species_init(378) = 0;
    species_init(379) = 0;
    species_init(380) = 0;
    species_init(381) = 0;
    species_init(382) = 0;
    species_init(383) = 0;
    species_init(384) = 0;
    species_init(385) = 0;
    species_init(386) = 0;
    species_init(387) = 0;
    species_init(388) = 0;
    species_init(389) = 0;
    species_init(390) = 0;
    species_init(391) = 0;
    species_init(392) = 0;
    species_init(393) = 0;
    species_init(394) = 0;
    species_init(395) = 0;
    species_init(396) = 0;
    species_init(397) = 0;
    species_init(398) = 0;
    species_init(399) = 0;
    species_init(400) = 0;
    species_init(401) = 0;
    species_init(402) = 0;
    species_init(403) = 0;
    species_init(404) = 0;
    species_init(405) = 0;
    species_init(406) = 0;
    species_init(407) = 0;
    species_init(408) = 0;
    species_init(409) = 0;
    species_init(410) = 0;
    species_init(411) = 0;
    species_init(412) = 0;
    species_init(413) = 0;
    species_init(414) = 0;
    species_init(415) = 0;
    species_init(416) = 0;
    species_init(417) = 0;
    species_init(418) = 0;
    species_init(419) = 0;
    species_init(420) = 0;
    species_init(421) = 0;
    species_init(422) = 0;
    species_init(423) = 0;
    species_init(424) = 0;
    species_init(425) = 0;
    species_init(426) = 0;
    species_init(427) = 0;
    species_init(428) = 0;
    species_init(429) = 0;
    species_init(430) = 0;
    species_init(431) = 0;
    species_init(432) = 0;
    species_init(433) = 0;
    species_init(434) = 0;
    species_init(435) = 0;
    species_init(436) = 0;
    species_init(437) = 0;
    species_init(438) = 0;
    species_init(439) = 0;
    species_init(440) = 0;
    species_init(441) = 0;
    species_init(442) = 0;
    species_init(443) = 0;
    species_init(444) = 0;
    species_init(445) = 0;
    species_init(446) = 0;
    species_init(447) = 0;
    species_init(448) = 0;
    species_init(449) = 0;
    species_init(450) = 0;
    species_init(451) = 0;
    species_init(452) = 0;
    species_init(453) = 0;
    species_init(454) = 0;
    species_init(455) = 0;
    species_init(456) = 0;
    species_init(457) = 0;
    species_init(458) = 0;
    species_init(459) = 0;
    species_init(460) = 0;
    species_init(461) = 0;
    species_init(462) = 0;
    species_init(463) = 0;
    species_init(464) = 0;
    species_init(465) = 0;
    species_init(466) = 0;
    species_init(467) = 0;
    species_init(468) = 0;
    species_init(469) = 0;
    species_init(470) = 0;
    species_init(471) = 0;
    species_init(472) = 0;
    species_init(473) = 0;
    species_init(474) = 0;
    species_init(475) = 0;
    species_init(476) = 0;
    species_init(477) = 0;
    species_init(478) = 0;
    species_init(479) = 0;
    species_init(480) = 0;
    species_init(481) = 0;
    species_init(482) = 0;
    species_init(483) = 0;
    species_init(484) = 0;
    species_init(485) = 0;
    species_init(486) = 0;
    species_init(487) = 0;
    species_init(488) = 0;
    species_init(489) = 0;
    species_init(490) = 0;
    species_init(491) = 0;
    species_init(492) = 0;
    species_init(493) = 0;
    species_init(494) = 0;
    species_init(495) = 0;
    species_init(496) = 0;
    species_init(497) = 0;
    species_init(498) = 0;
    species_init(499) = 0;
    species_init(500) = 0;
    species_init(501) = 0;
    species_init(502) = 0;
    species_init(503) = 0;
    species_init(504) = 0;
    species_init(505) = 0;
    species_init(506) = 0;
    species_init(507) = 0;
    species_init(508) = 0;
    species_init(509) = 0;
    species_init(510) = 0;
    species_init(511) = 0;
    species_init(512) = 0;
    species_init(513) = 0;
    species_init(514) = 0;
    species_init(515) = 0;
    species_init(516) = 0;
    species_init(517) = 0;
    species_init(518) = 0;
    species_init(519) = 0;
    species_init(520) = 0;
    species_init(521) = 0;
    species_init(522) = 0;
    species_init(523) = 0;
    species_init(524) = 0;
    species_init(525) = 0;
    species_init(526) = 0;
    species_init(527) = 0;
    species_init(528) = 0;
    species_init(529) = 0;
    species_init(530) = 0;
    species_init(531) = 0;
    species_init(532) = 0;
    species_init(533) = 0;
    species_init(534) = 0;
    species_init(535) = 0;
    species_init(536) = 0;
    species_init(537) = 0;
    species_init(538) = 0;
    species_init(539) = 0;
    species_init(540) = 0;
    species_init(541) = 0;
    species_init(542) = 0;
    species_init(543) = 0;
    species_init(544) = 0;
    species_init(545) = 0;
    species_init(546) = 0;
    species_init(547) = 0;
    species_init(548) = 0;
    species_init(549) = 0;
    species_init(550) = 0;
    species_init(551) = 0;
    species_init(552) = 0;
    species_init(553) = 0;
    species_init(554) = 0;
    species_init(555) = 0;
    species_init(556) = 0;
    species_init(557) = 0;
    species_init(558) = 0;
    species_init(559) = 0;
    species_init(560) = 0;
    species_init(561) = 0;
    species_init(562) = 0;
    species_init(563) = 0;
    species_init(564) = 0;
    species_init(565) = 0;
    species_init(566) = 0;
    species_init(567) = 0;
    species_init(568) = 0;
    species_init(569) = 0;
    species_init(570) = 0;
    species_init(571) = 0;
    species_init(572) = 0;
    species_init(573) = 0;
    species_init(574) = 0;
    species_init(575) = 0;
    species_init(576) = 0;
    species_init(577) = 0;
    species_init(578) = 0;
    species_init(579) = 0;
    species_init(580) = 0;
    species_init(581) = 0;
    species_init(582) = 0;
    species_init(583) = 0;
    species_init(584) = 0;
    species_init(585) = 0;
    species_init(586) = 0;
    species_init(587) = 0;
    species_init(588) = 0;
    species_init(589) = 0;
    species_init(590) = 0;
    species_init(591) = 0;
    species_init(592) = 0;
    species_init(593) = 0;
    species_init(594) = 0;
    species_init(595) = 0;
    species_init(596) = 0;
    species_init(597) = 0;
    species_init(598) = 0;
    species_init(599) = 0;
    species_init(600) = 0;
    species_init(601) = 0;
    species_init(602) = 0;
    species_init(603) = 0;
    species_init(604) = 0;
    species_init(605) = 0;
    species_init(606) = 0;
    species_init(607) = 0;
    species_init(608) = 0;
    species_init(609) = 0;
    species_init(610) = 0;
    species_init(611) = 0;
    species_init(612) = 0;
    species_init(613) = 0;
    species_init(614) = 0;
    species_init(615) = 0;
    species_init(616) = 0;
    species_init(617) = 0;
    species_init(618) = 0;
    species_init(619) = 0;
    species_init(620) = 0;
    species_init(621) = 0;
    species_init(622) = 0;
    species_init(623) = 0;
    species_init(624) = 0;
    species_init(625) = 0;
    species_init(626) = 0;
    species_init(627) = 0;
    species_init(628) = 0;
    species_init(629) = 0;
    species_init(630) = 0;
    species_init(631) = 0;
    species_init(632) = 0;
    species_init(633) = 0;
    species_init(634) = 0;
    species_init(635) = 0;
    species_init(636) = 0;
    species_init(637) = 0;
    species_init(638) = 0;
    species_init(639) = 0;
    species_init(640) = 0;
    species_init(641) = 0;
    species_init(642) = 0;
    species_init(643) = 0;
    species_init(644) = 0;
    species_init(645) = 0;
    species_init(646) = 0;
    species_init(647) = 0;
    species_init(648) = 0;
    species_init(649) = 0;
    species_init(650) = 0;
    species_init(651) = 0;
    species_init(652) = 0;
    species_init(653) = 0;
    species_init(654) = 0;
    species_init(655) = 0;
    species_init(656) = 0;
    species_init(657) = 0;
    species_init(658) = 0;
    species_init(659) = 0;
    species_init(660) = 0;
    species_init(661) = 0;
    species_init(662) = 0;
    species_init(663) = 0;
    species_init(664) = 0;
    species_init(665) = 0;
    species_init(666) = 0;
    species_init(667) = 0;
    species_init(668) = 0;
    species_init(669) = 0;
    species_init(670) = 0;
    species_init(671) = 0;
    species_init(672) = 0;
    species_init(673) = 0;
    species_init(674) = 0;
    species_init(675) = 0;
    species_init(676) = 0;
    species_init(677) = 0;
    species_init(678) = 0;
    species_init(679) = 0;
    species_init(680) = 0;
    species_init(681) = 0;
    species_init(682) = 0;
    species_init(683) = 0;
    species_init(684) = 0;
    species_init(685) = 0;
    species_init(686) = 0;
    species_init(687) = 0;
    species_init(688) = 0;
    species_init(689) = 0;
    species_init(690) = 0;
    species_init(691) = 0;
    species_init(692) = 0;
    species_init(693) = 0;
    species_init(694) = 0;
    species_init(695) = 0;
    species_init(696) = 0;
    species_init(697) = 0;
    species_init(698) = 0;
    species_init(699) = 0;
    species_init(700) = 0;
    species_init(701) = 0;
    species_init(702) = 0;
    species_init(703) = 0;
    species_init(704) = 0;
    species_init(705) = 0;
    species_init(706) = 0;
    species_init(707) = 0;
    species_init(708) = 0;
    species_init(709) = 0;
    species_init(710) = 0;
    species_init(711) = 0;
    species_init(712) = 0;
    species_init(713) = 0;
    species_init(714) = 0;
    species_init(715) = 0;
    species_init(716) = 0;
    species_init(717) = 0;
    species_init(718) = 0;
    species_init(719) = 0;
    species_init(720) = 0;
    species_init(721) = 0;
    species_init(722) = 0;
    species_init(723) = 0;
    species_init(724) = 0;
    species_init(725) = 0;

end


end
