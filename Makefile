#
#

FC = gfortran
# FC = f90
FFLAGS =

CWD = ./
MYL = ../lib/
MATH = dsp_math_MATRIX

OUT = dsp.exe

INC =
LIB =
# LIB =	-lmatmpp_sc -lMSL2

dsp:	$(CWD)dsp_header.f90\
	$(CWD)dsp_fbi_inhomo.f90\
	$(CWD)dsp_mag.f90\
	$(CWD)dsp_slv.f90\
	$(CWD)dsp_mi_couple.f90\
	$(CWD)dsp_main.f90
	
	$(FC) $(FFLAGS) -c $(CWD)dsp_header.f90
	$(FC) $(FFLAGS) -c $(CWD)dsp_fbi_inhomo.f90
	$(FC) $(FFLAGS) -c $(CWD)dsp_mag.f90
	$(FC) $(FFLAGS) -c $(CWD)dsp_slv.f90
	$(FC) $(FFLAGS) -c $(CWD)dsp_mi_couple.f90
	$(FC) $(FFLAGS) -c $(CWD)dsp_main.f90
	
	$(FC) $(FFLAGS) $(CWD)dsp_header.o\
			$(CWD)dsp_fbi_inhomo.o\
			$(CWD)dsp_mag.o\
			$(CWD)dsp_slv.o\
			$(CWD)dsp_mi_couple.o\
			$(CWD)dsp_main.o\
			-o $(OUT) $(LIB)


clean:
	rm -f ./*.o ./*.mod ./*.log ./*.exe ./*.err ./*.out

clear:
	rm -f ./*.o ./*.mod ./*.log ./*.err ./*.out

