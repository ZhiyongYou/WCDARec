#########################################
#  written by chen songzhan 2017/11/14  #
#                                       # 
#########################################
inc_dir = include/
src_dir = src/
obj_dir = obj/
bin_dir = bin/
tmp_dir = tmp/
CCOMPILER  = g++
DEBUG = -ggdb  #-g -gstabs -gstabs+ -ggdb -glevel

object = $(obj_dir)main.o \
         $(obj_dir)WCDAMcEvent.o $(obj_dir)WCDAMcRecEvent.o \
         $(obj_dir)LHCALEvent.o $(obj_dir)SaveEvent.o $(obj_dir)EventDict.o \
		 $(obj_dir)WCDARec.o $(obj_dir)CaliEvent.o

main: $(obj_dir)main.a
	$(CCOMPILER) -o main $(obj_dir)main.a `root-config --cflags --libs` -lMinuit -lMatrix

$(obj_dir)main.a: $(object)
	ar -r $(obj_dir)main.a $(object)

#main: $(object)
#	$(CCOMPILER) $(DEBUG) -o main $(object) `root-config --cflags --libs` -lMinuit -lMatrix -Wno-deprecated

$(obj_dir)main.o: main.cc Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)WCDAMcEvent.o: $(src_dir)WCDAMcEvent.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)WCDAMcRecEvent.o: $(src_dir)WCDAMcRecEvent.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)LHCALEvent.o: $(src_dir)LHCALEvent.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)SaveEvent.o: $(src_dir)SaveEvent.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)EventDict.o: $(src_dir)EventDict.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)WCDARec.o: $(src_dir)WCDARec.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

$(obj_dir)CaliEvent.o: $(src_dir)CaliEvent.cpp Makefile $(inc_dir)
	$(CCOMPILER) $(DEBUG) -c $< -o $@ -DSCAN -I $(inc_dir) `root-config --cflags --libs`

.PHONY : clean
 clean :
	rm main $(object)
