target remote | openocd -f board/stm32f429discovery.cfg -c "gdb_port pipe; log_output openocd.log"
monitor reset halt
load
