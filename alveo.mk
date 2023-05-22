
############################## Help Section ##############################
ifneq ($(findstring Makefile, $(MAKEFILE_LIST)), Makefile)
help:
	$(ECHO) "Makefile Usage:"
	$(ECHO) "  make all TARGET=<sw_emu/hw_emu/hw> PLATFORM=<FPGA platform>"
	$(ECHO) "      Command to generate the design for specified Target and Shell."
	$(ECHO) ""
	$(ECHO) "  make clean "
	$(ECHO) "      Command to remove the generated non-hardware files."
	$(ECHO) ""
	$(ECHO) "  make cleanall"
	$(ECHO) "      Command to remove all the generated files."
	$(ECHO) ""
	$(ECHO) "  make test PLATFORM=<FPGA platform>"
	$(ECHO) "      Command to run the application. This is same as 'run' target but does not have any makefile dependency."
	$(ECHO) ""
	$(ECHO) "  make run TARGET=<sw_emu/hw_emu/hw> PLATFORM=<FPGA platform>"
	$(ECHO) "      Command to run application in emulation."
	$(ECHO) ""
	$(ECHO) "  make build TARGET=<sw_emu/hw_emu/hw> PLATFORM=<FPGA platform>"
	$(ECHO) "      Command to build xclbin application."
	$(ECHO) ""
	$(ECHO) "  make host"
	$(ECHO) "      Command to build host application."
	$(ECHO) ""
endif

############################## Setting up Project Variables ##############################
TARGET := hw
include ./utils.mk

TEMP_DIR := ./_x.$(TARGET).$(XSA)
BUILD_DIR := ./build_dir.$(TARGET).$(XSA)

LINK_OUTPUT := $(BUILD_DIR)/ksw2_extz2_hw.link.xclbin
PACKAGE_OUT = ./package.$(TARGET)

VPP_PFLAGS := 
CMD_ARGS = $(BUILD_DIR)/ksw2_extz2_hw.xclbin
CXXFLAGS += -I$(XILINX_XRT)/include -I$(XILINX_VIVADO)/include -Wall -O0 -g -std=c++1y
LDFLAGS += -L$(XILINX_XRT)/lib -pthread -lOpenCL


########################## Checking if PLATFORM in allowlist #######################
PLATFORM_BLOCKLIST += nodma 
############################## Setting up Host Variables ##############################
#Include Required Host Source Files
CXXFLAGS += -Icommon/includes/xcl2
CXXFLAGS += -Icommon/includes/bitmap
CXXFLAGS += -Isubmodules/hls_arbitrary_precision_types/include  
CXXFLAGS += -Isubmodules/ksw2
HOST_SRCS += common/includes/xcl2/xcl2.cpp common/includes/bitmap/bitmap.cpp src/host.cpp submodules/ksw2/ksw2_extz2_sse.c
# Host compiler global settings
CXXFLAGS += -fmessage-length=0
LDFLAGS += -lrt -lstdc++ 

############################## Setting up Kernel Variables ##############################
# Kernel compiler global settings
VPP_FLAGS += -t $(TARGET) --platform $(PLATFORM) --save-temps

VPP_LDFLAGS_ksw2 = --config ./ksw2_extz2_hw.cfg
EXECUTABLE = ./host-ksw2
EMCONFIG_DIR = $(TEMP_DIR)

############################## Setting Targets ##############################
.PHONY: all clean cleanall docs emconfig
all: check-platform check-device check-vitis $(EXECUTABLE) $(BUILD_DIR)/ksw2_extz2_hw.xclbin emconfig

.PHONY: host
host: $(EXECUTABLE)

.PHONY: build
build: check-vitis check-device $(BUILD_DIR)/ksw2_extz2_hw.xclbin

.PHONY: xclbin
xclbin: build

############################## Setting Rules for Binary Containers (Building Kernels) ##############################
$(TEMP_DIR)/ksw2_extz2_hw.xo: src/ksw2_extz2_hw.cpp
	mkdir -p $(TEMP_DIR)
	v++ $(VPP_FLAGS) -c -k ksw2_extz2_hw --temp_dir $(TEMP_DIR)  -I'$(<D)' -o'$@' '$<'

$(BUILD_DIR)/ksw2_extz2_hw.xclbin: $(TEMP_DIR)/ksw2_extz2_hw.xo
	mkdir -p $(BUILD_DIR)
	v++ $(VPP_FLAGS) -l $(VPP_LDFLAGS) --temp_dir $(TEMP_DIR) $(VPP_LDFLAGS_ksw2) -o'$(LINK_OUTPUT)' $(+)
	v++ -p $(LINK_OUTPUT) $(VPP_FLAGS) --package.out_dir $(PACKAGE_OUT) -o $(BUILD_DIR)/ksw2_extz2_hw.xclbin

############################## Setting Rules for Host (Building Host Executable) ##############################
$(EXECUTABLE): $(HOST_SRCS) | check-xrt
	g++ -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

emconfig:$(EMCONFIG_DIR)/emconfig.json
$(EMCONFIG_DIR)/emconfig.json:
	emconfigutil --platform $(PLATFORM) --od $(EMCONFIG_DIR)

############################## Setting Essential Checks and Running Rules ##############################
run: all
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	cp -rf $(EMCONFIG_DIR)/emconfig.json .
	XCL_EMULATION_MODE=$(TARGET) $(EXECUTABLE) $(CMD_ARGS)
else
	$(EXECUTABLE) $(CMD_ARGS)
endif


.PHONY: test
test: $(EXECUTABLE)
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	XCL_EMULATION_MODE=$(TARGET) $(EXECUTABLE) $(CMD_ARGS)
else
	$(EXECUTABLE) $(CMD_ARGS)
endif


############################## Cleaning Rules ##############################
# Cleaning stuff
clean:
	-$(RMDIR) $(EXECUTABLE) $(XCLBIN)/{*sw_emu*,*hw_emu*} 
	-$(RMDIR) profile_* TempConfig system_estimate.xtxt *.rpt *.csv 
	-$(RMDIR) src/*.ll *v++* .Xil emconfig.json dltmp* xmltmp* *.log *.jou *.wcfg *.wdb

cleanall: clean
	-$(RMDIR) build_dir*
	-$(RMDIR) package.*
	-$(RMDIR) _x* *xclbin.run_summary qemu-memory-_* emulation _vimage pl* start_simulation.sh *.xclbin
