% fgraph-4topsSM.mp -- do not edit, generated automatically by hdr.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 110*1.0pt#, 60*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleftn (__i, 2);
vrightn (__o, 4);
save curly_len; curly_len:=2mm;
vlabel ("$g$", __i1);
vlabel ("$g$", __i2);
vlabel ("$t$", __o1);
vlabel ("$\overline {t}$", __o2);
vlabel ("$t$", __o3);
vlabel ("$\overline {t}$", __o4);
vconnect ("gluon", __i1, __v1);
vconnect ("fermion", __v1, __o1);
vconnect ("fermion", __v2, __v1);
vconnect ("fermion", __o2, __v4);
vconnect ("gluon", __i2, __v3);
vconnect ("fermion", __o4, __v3);
vconnect ("fermion", __v3, __v2);
vconnect ("fermion", __v4, __o3);
vconnect ("gluon", __v2, __v4);
vconstraint (((40),whatever), __i1, __v1);
vconstraint (((40),whatever), __i2, __v3);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
