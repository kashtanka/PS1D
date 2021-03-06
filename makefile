#makefile for the single-column model
model_path= /home/dmitry/models/PS_1D/
source_path = /home/dmitry/models/PS_1D/source
obj_path = /home/dmitry/models/PS_1D/obj
exec = ps_1d.out
fc = ifort
switch = -r8 -O3 -traceback
obj = \
 $(obj_path)/alloc_1d_mod.o \
 $(obj_path)/ice_mod.o \
 $(obj_path)/goddard_mod.o \
 $(obj_path)/alloc_goddard.o \
 $(obj_path)/cliradlw.o \
 $(obj_path)/ice.o \
 $(obj_path)/ps1d_main.o \
 $(obj_path)/readpa.o \
 $(obj_path)/vertical_grid.o \
 $(obj_path)/allocvar.o \
 $(obj_path)/iniprofs.o \
 $(obj_path)/ini_ice.o \
 $(obj_path)/eq_state.o \
 $(obj_path)/radiation.o \
 $(obj_path)/moment.o \
 $(obj_path)/vert_adv.o \
 $(obj_path)/thermo.o \
 $(obj_path)/humid.o \
 $(obj_path)/bl_depth.o \
 $(obj_path)/aselyn.o \
 $(obj_path)/baroclinicity.o \
 $(obj_path)/surf_layer_t.o \
 $(obj_path)/balance.o \
 $(obj_path)/const_fluxes.o \
 $(obj_path)/diffu_local.o \
 $(obj_path)/richardson.o \
 $(obj_path)/diffu_LS96.o \
 $(obj_path)/diffu_Noh03.o \
 $(obj_path)/diffu_TM86.o \
 $(obj_path)/diffu_INM.o \
 $(obj_path)/diffu_Lock.o \
 $(obj_path)/diffu_ml.o \
 $(obj_path)/implicit_dif.o \
 $(obj_path)/progonka.o \
 $(obj_path)/microphysics.o \
 $(obj_path)/microphys_update.o \
 $(obj_path)/satur_adjust.o \
 $(obj_path)/entrainment.o \
 $(obj_path)/crossection.o \
 $(obj_path)/meanABL.o \
 $(obj_path)/output.o \
#
$(exec) : $(obj)
	$(fc) -o $(exec) $(switch) $(obj)
#MODULES:
$(obj_path)/alloc_1d_mod.o : $(source_path)/alloc_1d_mod.for
	cd $(obj_path)/ && $(fc) -c $(switch) $(source_path)/alloc_1d_mod.for && cd $(model_path)
$(obj_path)/ice_mod.o : $(source_path)/ice_mod.f
	cd $(obj_path)/ && $(fc) -c $(switch) $(source_path)/ice_mod.f && cd $(model_path)
#OTHER SOURCE
$(obj_path)/%.o : $(source_path)/%.for
	cd $(obj_path)/ && $(fc) -c $(switch) $< && cd $(model_path)
$(obj_path)/%.o : $(source_path)/%.f
	cd $(obj_path)/ && $(fc) -c $(switch) $< && cd $(model_path)
clean :
	rm $(obj) $(exec)

