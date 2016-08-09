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
tmpfs_size = $(shell echo "$$(du -L --max-depth 0 $(read_folder) | sed -r "s/\t.+//")*2000" | bc )

.PHONY: all trimgalore trimgalore_stats 2_filter_pe 2_filter_unpaired 2_filter_stats

all: 2_filter_stats

init: 
	mkdir -p tmp
	sudo mount -t tmpfs -o size=$(tmpfs_size) tmpfs $(tmpfs_mount_folder)

trimgalore: init 
	mkdir -p log
	time 1_run_trimgalore.sh $(read_folder) $(tmpfs_mount_folder) 2>&1 | tee log/${sample_id}_$@.step.log

trimgalore_stats: trimgalore
	zcat $(tmpfs_mount_folder)/*_val_1.fq.gz | sed -n "2~4p" | awk 'BEGIN{bases=0;reads=0;} {bases+=length($1);reads+=1} END{print reads,bases}' > log/${sample_id}_trimgalore_R1.stats
	zcat $(tmpfs_mount_folder)/*_val_2.fq.gz | sed -n "2~4p" | awk 'BEGIN{bases=0;reads=0;} {bases+=length($1);reads+=1} END{print reads,bases}' > log/${sample_id}_trimgalore_R2.stats
	zcat $(tmpfs_mount_folder)/*_unpaired_?.fq.gz | sed -n "2~4p" | awk 'BEGIN{bases=0;reads=0;} {bases+=length($1);reads+=1} END{print reads,bases}' > log/${sample_id}_trimgalore_unpaired.stats

#2nd step: Remove SISPA overhangs and discard human reads
2_filter_pe: trimgalore
	mkdir -p log
	time 2_filter_pe.sh $(sample_id) $(tmpfs_mount_folder) 2>&1 | tee log/${sample_id}_$@.step.log

2_filter_unpaired: trimgalore_stats
	mkdir -p log
	time 2_filter_unpaired.sh $(sample_id) $(tmpfs_mount_folder) 2>&1 | tee log/${sample_id}_$@.step.log

2_filter_stats: 2_filter_pe 2_filter_unpaired
	zcat $(sample_id)_pe.fq.gz | sed -n "2~4p" | awk 'BEGIN{bases=0;reads=0;} {bases+=length($1);reads+=1} END{print reads,bases}' > log/${sample_id}_2_filter_pe.stats
	zcat $(sample_id)_unpaired.fq.gz | sed -n "2~4p" | awk 'BEGIN{bases=0;reads=0;} {bases+=length($1);reads+=1} END{print reads,bases}' > log/${sample_id}_2_filter_unpaired.stats

clean:
	rm -i $(tmpfs_mount_folder)/*
	sudo umount $(tmpfs_mount_folder)
