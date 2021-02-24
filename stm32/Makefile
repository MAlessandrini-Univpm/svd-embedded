# Makefile to build firmware with gcc on linux
# author: Michele Alessandrini

KALMAN = 0

BASE_DIR = .
COMMON_DIR = ./f4xx_common
PROGRAMNAME = svd_stm32
RELEASE ?= 1
ifeq ($(RELEASE), 1)
BUILD_DIR = ./Release
else
BUILD_DIR = ./Debug
endif

GCC_EXTRA_OPTIONS = \
	-DUSE_STDPERIPH_DRIVER -DSTM32F429_439xx -DUSE_STM32F429I_DISCO \
	-I$(BASE_DIR) \
	-I$(COMMON_DIR)/Libraries/STM32F4xx_StdPeriph_Driver/inc \
	-I$(COMMON_DIR)/Libraries/CMSIS/Include \
	-I$(COMMON_DIR)/Libraries/CMSIS/Device/ST/STM32F4xx/Include \
	#-I$(COMMON_DIR)/Utilities/STM32F429I-Discovery

LIBRARIES =

SRCS_C_BASE = $(wildcard $(BASE_DIR)/*.c)
ifeq ($(KALMAN), 1)
SRCS_C_BASE := $(filter-out $(BASE_DIR)/main.c, $(SRCS_C_BASE))
else
SRCS_C_BASE := $(filter-out $(BASE_DIR)/main_kalman.c $(BASE_DIR)/kalman.c, $(SRCS_C_BASE))
endif

SRCS_C_COMMON = \
	$(COMMON_DIR)/Libraries/STM32F4xx_StdPeriph_Driver/src/stm32f4xx_gpio.c \
	$(COMMON_DIR)/Libraries/STM32F4xx_StdPeriph_Driver/src/stm32f4xx_rcc.c \
#	$(COMMON_DIR)/Libraries/STM32F4xx_StdPeriph_Driver/src/stm32f4xx_dma2d.c \
#	$(COMMON_DIR)/Libraries/STM32F4xx_StdPeriph_Driver/src/stm32f4xx_spi.c \
#	$(COMMON_DIR)/Libraries/STM32F4xx_StdPeriph_Driver/src/stm32f4xx_exti.c \
#	$(COMMON_DIR)/Libraries/STM32F4xx_StdPeriph_Driver/src/stm32f4xx_i2c.c \
#	$(COMMON_DIR)/Libraries/STM32F4xx_StdPeriph_Driver/src/stm32f4xx_ltdc.c \
#	$(COMMON_DIR)/Libraries/STM32F4xx_StdPeriph_Driver/src/stm32f4xx_fmc.c \
#	$(COMMON_DIR)/Utilities/STM32F429I-Discovery/stm32f429i_discovery_sdram.c \
#	$(COMMON_DIR)/Utilities/STM32F429I-Discovery/stm32f429i_discovery_ioe.c \
#	$(COMMON_DIR)/Utilities/STM32F429I-Discovery/stm32f429i_discovery_lcd.c

SRCS_ASM_BASE = $(wildcard $(BASE_DIR)/*.s)
SRCS_ASM_COMMON = $(COMMON_DIR)/startup_stm32f429_439xx.s

OBJS  = $(addprefix $(BUILD_DIR)/, $(patsubst $(BASE_DIR)/%, %, $(addsuffix .o, $(basename $(SRCS_C_BASE) $(SRCS_ASM_BASE)))))
OBJS += $(addprefix $(BUILD_DIR)/, $(patsubst $(COMMON_DIR)/%, %, $(addsuffix .o, $(basename $(SRCS_C_COMMON) $(SRCS_ASM_COMMON)))))
DEPS  = $(addprefix $(BUILD_DIR)/, $(patsubst $(BASE_DIR)/%, %, $(addsuffix .d, $(basename $(SRCS_C_BASE)))))
DEPS += $(addprefix $(BUILD_DIR)/, $(patsubst $(COMMON_DIR)/%, %, $(addsuffix .d, $(basename $(SRCS_C_COMMON)))))


COMPILER = arm-none-eabi-gcc
LINKER = arm-none-eabi-ld
OBJCOPY = arm-none-eabi-objcopy
OBJDUMP = arm-none-eabi-objdump
GDB = gdb-multiarch

CXXFLAGS = -mcpu=cortex-m4 -mthumb -mlittle-endian -mfloat-abi=hard -mfpu=fpv4-sp-d16 \
	-Wall -Wextra -Wno-missing-field-initializers -std=c99 -g

ifeq ($(RELEASE), 1)
CXXFLAGS += -O3
endif

CXXFLAGS += $(GCC_EXTRA_OPTIONS)

LINKSCRIPT = stm32f429zi.ld
GCCLINKFLAGS = -static -Wl,--entry=Reset_Handler --specs=nosys.specs

#ifeq ($(RELEASE), 1)
#GCCLINKFLAGS += -Wl,--strip-debug
#endif


TARGETS = $(BUILD_DIR)/$(PROGRAMNAME).out $(BUILD_DIR)/$(PROGRAMNAME).srec $(BUILD_DIR)/$(PROGRAMNAME).bin $(BUILD_DIR)/$(PROGRAMNAME).out.dump


# targets

all: $(TARGETS)

clean:
	rm -f $(TARGETS) $(OBJS) $(DEPS) $(BUILD_DIR)/*.map
	find $(BUILD_DIR) -depth -type d -empty -exec rmdir {} \;

cleanTargets:
	rm -f $(TARGETS)

# rules

include $(DEPS)

$(BUILD_DIR)/%.o: $(BASE_DIR)/%.c
	-mkdir -p $(dir $@)
	$(COMPILER) $(CXXFLAGS) -o $@ -c $<

$(BUILD_DIR)/%.o: $(COMMON_DIR)/%.c
	-mkdir -p $(dir $@)
	$(COMPILER) $(CXXFLAGS) -o $@ -c $<

$(BUILD_DIR)/%.o: $(BASE_DIR)/%.S
	-mkdir -p $(dir $@)
	$(COMPILER) $(CXXFLAGS) -o $@ -c $<

$(BUILD_DIR)/%.o: $(COMMON_DIR)/%.S
	-mkdir -p $(dir $@)
	$(COMPILER) $(CXXFLAGS) -o $@ -c $<

$(BUILD_DIR)/%.o: $(BASE_DIR)/%.s
	-mkdir -p $(dir $@)
	$(COMPILER) $(CXXFLAGS) -o $@ -c $<

$(BUILD_DIR)/%.o: $(COMMON_DIR)/%.s
	-mkdir -p $(dir $@)
	$(COMPILER) $(CXXFLAGS) -o $@ -c $<


# automatic C dependency handling
$(BUILD_DIR)/%.d: $(BASE_DIR)/%.c
	-mkdir -p $(dir $@)
	$(COMPILER) -MM -MT $(addsuffix .o, $(basename $@)) $(CXXFLAGS) $< | sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $@

$(BUILD_DIR)/%.d: $(COMMON_DIR)/%.c
	-mkdir -p $(dir $@)
	$(COMPILER) -MM -MT $(addsuffix .o, $(basename $@)) $(CXXFLAGS) $< | sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $@


%.srec: %.out
	$(OBJCOPY) -O srec  $< $@

%.bin: %.out
	$(OBJCOPY) -O binary  $< $@

%.out.dump: %.out
	$(OBJDUMP) -dxsSt $< > $@



$(BUILD_DIR)/%.out: $(OBJS)
	$(COMPILER) -T $(LINKSCRIPT) $(CXXFLAGS) $(GCCLINKFLAGS) -o $@ $^ $(LIBRARIES) -Wl,-Map=$@.map


flash: $(BUILD_DIR)/$(PROGRAMNAME).out
	openocd -f board/stm32f429discovery.cfg -c "program $< verify reset exit"


dbg: $(BUILD_DIR)/$(PROGRAMNAME).out .gdbinit
	$(GDB) -tui $(BUILD_DIR)/$(PROGRAMNAME).out


# target to generate files for Qt Creator IDE
IDE_FILES = $(SRCS_C_BASE) $(SRCS_C_COMMON) $(SRCS_ASM_BASE) $(SRCS_ASM_COMMON) $(wildcard $(BASE_DIR)/*.h) Makefile
qtcreator_files:
	@echo $(IDE_FILES) | tr ' ' '\n' | sort > $(PROGRAMNAME).files
	@echo $(CXXFLAGS) | tr ' ' '\n' | grep '^-I' | sed 's/^-I//' > $(PROGRAMNAME).includes
	@echo $(CXXFLAGS) | tr ' ' '\n' | grep '^-D' | sed 's/^-D/#define /' | tr '=' ' ' > $(PROGRAMNAME).config
	@touch $(PROGRAMNAME).creator


.PHONY: all clean flash dbg
.SECONDARY: $(OBJS)
FORCE:

