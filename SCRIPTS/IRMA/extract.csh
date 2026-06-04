foreach n(00)
echo $n
#ncrcat -v USFC,VSFC,PSFC,PW,TB -d lon,1700,1950 -d lat,630,880 ../../OUT_2D/KATRINA_2304x1152x74_15s_mom23_0.8_ens{$n}_0.4K_1152_20*0.2D_atm.nc ../../OUT_2D/KATRINA_2304x1152x74_15s_mom23_0.8_ens{$n}_0.4K_KATRINA.2D_atm.nc
ncrcat -O -v  USFC,VSFC,PSFC,PW,TB -d lon,1750,2100 -d lat,630,880 ../../OUT_2D/IRMA_2304x1152x74_15s_Sep5_mom23_dampw5u5_fudge0.3_new_{$n}_1152_201709*0.2D_atm.nc ../../OUT_2D/IRMA_2304x1152x74_15s_Sep5_mom23_dampw5u5_fudge0.3_new_{$n}_IRMA.2D_atm.nc
end
