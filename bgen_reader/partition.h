void read_partition_part1(struct bgen_partition const* partition, uint32_t* position,
                          uint16_t* nalleles, uint64_t* offset, uint32_t* id_max_len,
                          uint32_t* rsid_max_len, uint32_t* chrom_max_len,
                          uint32_t* allele_ids_max_len);

void read_partition_part2(struct bgen_partition const* partition, wchar_t* const id,
                          uint32_t id_stride, wchar_t* const rsid, uint32_t rsid_stride,
                          wchar_t* const chrom, uint32_t chrom_stride, char* const allele_ids,
                          uint32_t allele_ids_stride);
