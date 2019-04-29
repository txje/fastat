#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include "incl/klib/kseq.h"
#include "incl/klib/kvec.h"
#include "incl/klib/ksort.h"
#include <getopt.h>
#include <stdint.h>

//#include "incl/svg.h"

KSORT_INIT_GENERIC(int)

KSEQ_INIT(gzFile, gzread)

void usage() {
  printf("Usage: fastat [options] FASTA/Q[.gz] ...\n");
  printf("Options:\n");
  printf("  -v, --verbose: verbose\n");
  printf("  -j, --joint: do aggregate stats over all files\n");
  printf("  -c, --compact: print parseable compact stats\n");
  printf("  -a, --all print size/GC stats for all sequences\n");
  printf("  -h, --help: show this\n");
}

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "verbose",                no_argument,       0, 'v' },
  { "joint",                  no_argument,       0, 'j' },
  { "compact",                no_argument,       0, 'c' },
  { "all",                    no_argument,       0, 'a' },
  { "help",                   no_argument,       0, 'h' },
  { 0, 0, 0, 0}
};

/*
int plot() {
  svg* psvg;
  psvg = svg_create(192, 320);

  if(psvg == NULL) {
    puts("psvg is NULL");
  } else {
    svg_fill(psvg, "#DADAFF");

    svg_text(psvg, 32, 32, "sans-serif", 16, "#000000", "#000000", "drawallshapes");
    svg_circle(psvg, "#000080", 4, "#0000FF", 32, 64, 96);
    svg_ellipse(psvg, 64, 160, 32, 16, "#FF0000", "#800000", 4);

    svg_line(psvg, "#000000", 2, 32, 192, 160, 192);

    svg_rectangle(psvg, 64, 64, 32, 224, "#00FF00", "#008000", 4, 4, 4);

    svg_finalize(psvg);
    svg_print(psvg);
    svg_save(psvg, "allshapes.svg");
    svg_free(psvg);
  }
}
*/

int main(int argc, char *argv[]) {
  char* fasta = NULL;
  int verbose = 0;
  int joint = 0;
  int compact = 0;
  int all = 0;

  int opt, long_idx, index;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "vhjca", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'j':
        joint = 1;
        break;
      case 'c':
        compact = 1;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'a':
        all = 1;
        break;
      case 'h':
        usage();
        return 0;
        break;
      default:
        usage();
        return 1;
    }
  }

  gzFile f;
  kseq_t *ks;

  uint64_t an = 0;
  uint64_t atotal_size = 0;
  kvec_t(int) aa;
  kv_init(aa);

  int i; // a reusable counter

  if(optind >= argc) { // no positional arguments (FASTAs)
    fprintf(stderr, "FASTA/Q[.gz] is required\n");
    usage();
    return 1;
  }

  if(compact)
    fprintf(stdout, "file\tn_seqs\ttotal_size\tavg_size\tmedian\tmaximum\tN50\n");

  for (index = optind; index < argc; index++) {
    fasta = argv[index];

    f = gzopen(fasta, "r");
    if(!f) {
      fprintf(stderr, "'%s' does not exist, try again\n\n", fasta);
      usage();
      return 1;
    }
    ks = kseq_init(f); 

    uint64_t n = 0;
    uint64_t total_size = 0;
    uint64_t gc = 0;
    kvec_t(int) a;
    kv_init(a);
    while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
      if(verbose) {
        fprintf(stderr, "Processing read %d (%s, %u bp): %s\n", n, ks->name.s, ks->seq.l, ks->seq.s);
      }
      total_size += ks->seq.l;
      if(all) {
        gc = 0;
        for(i = 0; i < ks->seq.l; i++)
          if(ks->seq.s[i] == 'G' || ks->seq.s[i] == 'C' || ks->seq.s[i] == 'g' || ks->seq.s[i] == 'c')
            gc++;
        fprintf(stdout, "%s\t%u\t%f\n", ks->name.s, ks->seq.l, ((float)gc)/ks->seq.l);
      }
      n++;
      kv_push(int, a, ks->seq.l);

      atotal_size += ks->seq.l;
      an++;
      if(joint) {
        kv_push(int, aa, ks->seq.l);
      }
    }
    ks_mergesort(int, kv_size(a), a.a, 0);

    if(!joint) {
      int i;

      if(!compact) {
        fprintf(stdout, "%s\n", fasta);
        fprintf(stdout, "%llu reads, %llu bp, mean: %f, median: %d\n", n, total_size, (double)total_size/n, kv_A(a, kv_size(a)/2));

        if(kv_size(a) > 0) {
          fprintf(stdout, "Largest:");
          for(i = 1; i <= (5 >= kv_size(a) ? kv_size(a) : 5); i++) {
            if(i > 1) fprintf(stdout, ",");
            fprintf(stdout, " %u", kv_A(a, kv_size(a)-i));
          }
          fprintf(stdout, "\n");
        }
      } else {
        fprintf(stdout, "%s\t%llu\t%llu\t%f\t%d\t%d\t", fasta, n, total_size, (double)total_size/n, kv_A(a, kv_size(a)/2), kv_A(a, kv_size(a)-1));
      }

      uint64_t cum = 0;
      for(i = kv_size(a)-1; i >= 0; i--) {
        if((cum * 10) / total_size != ((cum + (uint32_t)kv_A(a, i)) * 10) / total_size) {
          if(!compact) {
            fprintf(stdout, "N%llu: %d\n", ((cum + (uint32_t)kv_A(a, i)) * 10 / total_size) * 10, kv_A(a, i));
          }
          else if(cum*10/total_size < 5 && (cum + (uint32_t)kv_A(a,i))*10/total_size >= 5)
            fprintf(stdout, "%d\n", kv_A(a, i));
        }
        cum += kv_A(a, i);
      }
    }

    kv_destroy(a);
    kseq_destroy(ks);
    gzclose(f);
  }

  if(joint) {
    ks_mergesort(int, kv_size(aa), aa.a, 0);
    int i;

    if(!compact) {
      fprintf(stdout, "Total:\n");
      fprintf(stdout, "%llu reads, %llu bp, mean: %f, median: %d\n", an, atotal_size, (double)atotal_size/an, kv_A(aa, kv_size(aa)/2));

      if(kv_size(aa) > 0) {
        fprintf(stdout, "Largest:");
        for(i = 1; i <= (5 >= kv_size(aa) ? kv_size(aa) : 5); i++) {
          if(i > 1) fprintf(stdout, ",");
          fprintf(stdout, " %u", kv_A(aa, kv_size(aa)-i));
        }
        fprintf(stdout, "\n");
      }
    } else {
      fprintf(stdout, "%s\t%llu\t%llu\t%f\t%d\t%d\t", "total", an, atotal_size, (double)atotal_size/an, kv_A(aa, kv_size(aa)/2), kv_A(aa, kv_size(aa)-1));
    }

    uint64_t cum = 0;
    for(i = kv_size(aa)-1; i >= 0; i--) {
      if((cum * 10) / atotal_size != ((cum + (uint32_t)kv_A(aa, i)) * 10) / atotal_size) {
        if(!compact) {
          fprintf(stdout, "N%llu: %d\n", ((cum + (uint32_t)kv_A(aa, i)) * 10 / atotal_size) * 10, kv_A(aa, i));
        }
        else if(cum*10/atotal_size < 5 && (cum + (uint32_t)kv_A(aa,i))*10/atotal_size >= 5)
          fprintf(stdout, "%d\n", kv_A(aa, i));
      }
      cum += kv_A(aa, i);
    }
  }

  return 0;
}
