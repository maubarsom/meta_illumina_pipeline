# Makefile for extracting filtered microbiome reads

ifndef sample_id
$(error Define sample_id variable)
endif

ifndef read_folder
read_folder := reads
$(warning 'read_folder' variable defaults to $(read_folder))
endif

tmpfs_mount_folder := ./tmp
#Twice the size of the raw reads
tmpfs_size = $(shell echo "$$(du -L $(read_folder) | sed -r "s/\t.+//")*2000" | bc )

.PHONY: all trimgalore 2_filter_pe

all: 2_filter_pe

init: 
	mkdir -p tmp
	sudo mount -t tmpfs -o size=$(tmpfs_size) tmpfs $(tmpfs_mount_folder)

trimgalore: init 
	time bash 1_run_trimgalore.sh $(read_folder) $(tmpfs_mount_folder) 2>&1 | tee trimgalore.step.log

#trimgalore_stats: trimgalore
#	time bash 2_trimgalore_stats.sh $(tmpfs_mount_folder)

#2nd step: Remove SISPA overhangs and discard human reads
2_filter_pe: trimgalore
	time bash 2_filter_pe.sh $(sample_id) $(tmpfs_mount_folder) 

#2_filter_se: trimgalore
#	time bash 2_filter_se.sh $(sample_id) $(tmpfs_mount_folder)

clean:
	rm -i $(tmpfs_mount_folder)/*
	sudo umount $(tmpfs_mount_folder)
