CC = gcc
CFLAGS = -fPIC -std=gnu99 -O3 -march=native -flto -minline-stringops-dynamically -g
CPPFLAGS = -Iinchi -DTARGET_API_LIB -I$(JAVA_HOMEX)/include -I$(JAVA_HOMEX)/include/linux
LDFLAGS = -shared
TARGET_LIB = libelchem.so

SRCS = \
	inchi/ichi_bns.c \
	inchi/ichicano.c \
	inchi/ichicans.c \
	inchi/ichican2.c \
	inchi/ichiisot.c \
	inchi/ichimake.c \
	inchi/ichimap1.c \
	inchi/ichimap2.c \
	inchi/ichimap4.c \
	inchi/ichinorm.c \
	inchi/ichiqueu.c \
	inchi/ichiring.c \
	inchi/ichisort.c \
	inchi/ichister.c \
	inchi/ichitaut.c \
	inchi/readinch.c \
	inchi/strutil.c \
	inchi/util.c \
	inchi.c \
	isomorphism.c

OUTDIR ?= .
DEPDIR = .deps
OBJDIR = .objs

OBJS := $(patsubst %,$(OBJDIR)/%.o,$(notdir $(basename $(SRCS))))
DEPS := $(patsubst %,$(DEPDIR)/%.d,$(notdir $(basename $(SRCS))))

DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
COMPILE.c = $(CC) $(DEPFLAGS) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c

.PHONY: all claen

all: $(OUTDIR)/$(TARGET_LIB)

$(OUTDIR)/$(TARGET_LIB): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJDIR)/%.o : %.c $(DEPDIR)/%.d
	$(COMPILE.c) $(OUTPUT_OPTION) $<

$(OBJDIR)/%.o : inchi/%.c $(DEPDIR)/%.d
	$(COMPILE.c) $(OUTPUT_OPTION) $<

$(DEPDIR):
	@mkdir -p $@

$(OBJDIR):
	@mkdir -p $@

clean:
	rm -f $(DEPDIR)/*.d $(OBJDIR)/*.o $(OUTDIR)/$(TARGET_LIB) 
	if [ -d "$(OBJDIR)" ]; then rmdir $(OBJDIR); fi
	if [ -d "$(DEPDIR)" ]; then rmdir $(DEPDIR); fi

$(OBJS): | $(OBJDIR)

$(DEPS): | $(DEPDIR) 

$(DEPS):
include $(DEPS)
