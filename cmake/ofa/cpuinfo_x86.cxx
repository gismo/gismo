#include <cstdio>
#include <cstring>
#include <string>

#define print_features(reg,features,n)                      \
  for (int i=0; i<n; ++i)                                   \
    printf ("%s", (reg>>i & 0x1) && !features[i].empty()    \
            ? (features[i]+" ").c_str() : "");

// Get the vendor ID
void getVendorID() {
  int32_t a[3];
  for(int i=0; i<3; ++i)
    a[i] = 0;
  
  // EAX=0: Vendor ID
  __asm__("mov $0x0, %eax\n\t");
  __asm__("cpuid\n\t");
  __asm__("mov %%ebx, %0\n\t":"=r" (a[0]));
  __asm__("mov %%edx, %0\n\t":"=r" (a[1]));
  __asm__("mov %%ecx, %0\n\t":"=r" (a[2]));

  char vendorID[13]; vendorID[12] = 0;
  memcpy(&vendorID[0],&a[0],4);
  memcpy(&vendorID[4],&a[1],4);
  memcpy(&vendorID[8],&a[2],4);
  
  printf ("vendor_id       : %s\n", vendorID);
}

// Get processor information
void getProcInfo() {
  int32_t eax = 0;

  // EAX=1: Processor Info
  __asm__("mov $0x1 , %eax\n\t");
  __asm__("cpuid\n\t");
  __asm__("mov %%eax, %0\n\t":"=r" (eax)); //gives model and family

  int32_t stepping = eax>>0 & 0xF;
  int32_t model    = eax>>4 & 0xF;
  int32_t family   = eax>>8 & 0xF;
  if(family == 6 || family == 15)
    model += (eax>>16 & 0xF)<<4;

  printf ("cpu family      : %d\n", family);
  printf ("model           : %d\n", model);
  printf ("stepping        : %d\n", stepping);
}

// Get processor features
void getFeatures() {
  int32_t eax_max,ecx_max,eax,ebx,ecx,edx;

  // CPU flags
  printf ("flags           : ");
  
  // EAX=0: largest value that EAX can be set to before calling CPUID
  __asm__("mov $0x0, %eax\n\t");
  __asm__("cpuid\n\t");
  __asm__("mov %%eax, %0\n\t":"=r" (eax_max));

  printf("\neax_max:%d\n", eax_max);
  
  if (eax_max >= 1) {
      
    // EAX=1: Processor Info and Feature Bits
    __asm__("mov $0x1 , %eax\n\t");
    __asm__("cpuid\n\t");
    __asm__("mov %%ecx, %0\n\t":"=r" (ecx)); //feature flags
    __asm__("mov %%edx, %0\n\t":"=r" (edx)); //feature flags   
      
    // Intel-defined CPU features, CPUID level 0x00000001 (EDX), word 0
    {
      std::string features[] = { "fpu",     /* Onboard FPU */
                                 "vme",     /* Virtual Mode Extensions */
                                 "de",      /* Debugging Extensions */
                                 "pse",     /* Page Size Extensions */
                                 "tsc",     /* Time Stamp Counter */
                                 "msr",     /* Model-Specific Registers */
                                 "pae",     /* Physical Address Extensions */
                                 "mce",     /* Machine Check Exception */
                                 "cx8",     /* CMPXCHG8 instruction */
                                 "apic",    /* Onboard APIC */
                                 "",        /* Reserved */
                                 "sep",     /* SYSENTER/SYSEXIT */
                                 "mtrr",    /* Memory Type Range Registers */
                                 "pge",     /* Page Global Enable */
                                 "mca",     /* Machine Check Architecture */
                                 "cmov",    /* CMOV instructions (plus FCMOVcc, FCOMI with FPU) */
                                 "pat",     /* Page Attribute Table */
                                 "pse36",   /* 36-bit PSEs */
                                 "pn",      /* Processor serial number */
                                 "clflush", /* CLFLUSH instruction */
                                 "",        /* Reserved */
                                 "dts",     /* "dts" Debug Store */
                                 "acpi",    /* ACPI via MSR */
                                 "mmx",     /* Multimedia Extensions */
                                 "fxsr",    /* FXSAVE/FXRSTOR, CR4.OSFXSR */
                                 "sse",     /* "sse" */
                                 "sse2",    /* "sse2" */
                                 "ss",      /* "ss" CPU self snoop */
                                 "ht",      /* Hyper-Threading */
                                 "tm",      /* "tm" Automatic clock control */
                                 "ia64",    /* IA-64 processor */
                                 "pbe"      /* Pending Break Enable */
      };
      printf("INTEL 0x00000001 (EDX)");
      print_features(edx, features, 32);
    }
    
    // Intel-defined CPU features, CPUID level 0x00000001 (ECX), word 4
    {
      std::string features[] = { "sse3",      /* "pni" SSE-3 */
                                 "pclmulqdq", /* PCLMULQDQ instruction */
                                 "dtes64",    /* 64-bit Debug Store */
                                 "monitor",   /* "monitor" MONITOR/MWAIT support */
                                 "ds_cpl",    /* "ds_cpl" CPL-qualified (filtered) Debug Store */
                                 "vmx",       /* Hardware virtualization */
                                 "smx",       /* Safer Mode eXtensions */
                                 "est",       /* Enhanced SpeedStep */
                                 "tm2",       /* Thermal Monitor 2 */
                                 "ssse3",     /* Supplemental SSE-3 */
                                 "cid",       /* Context ID */
                                 "sdbg",      /* Silicon Debug */
                                 "fma",       /* Fused multiply-add */
                                 "cx16",      /* CMPXCHG16B instruction */
                                 "xtpr",      /* Send Task Priority Messages */
                                 "pdcm",      /* Perf/Debug Capabilities MSR */
                                 "",          /* Reserved */
                                 "pcid",      /* Process Context Identifiers */
                                 "dca",       /* Direct Cache Access */
                                 "sse4_1",    /* "sse4_1" SSE-4.1 */
                                 "sse4_2",    /* "sse4_2" SSE-4.2 */
                                 "x2apic",    /* X2APIC */
                                 "movbe",     /* MOVBE instruction */
                                 "popcnt",    /* POPCNT instruction */
                                 "tsc_deadline_timer", /* TSC deadline timer */
                                 "aes",       /* AES instructions */
                                 "xsave",     /* XSAVE/XRSTOR/XSETBV/XGETBV instructions */
                                 "",          /* "" XSAVE instruction enabled in the OS */
                                 "avx",       /* Advanced Vector Extensions */
                                 "f16c",      /* 16-bit FP conversions */
                                 "rdrand",    /* RDRAND instruction */
                                 "hypervisor" /* Running on a hypervisor */
      };
      printf("INTEL 0x00000001 (ECX)");
      print_features(ecx, features, 32);
    }
  } // EAX=1
  
  if (eax_max >= 7) {      
    // EAX=7, ECX=0: Extended Features
    __asm__("mov $0x7 , %eax\n\t");
    __asm__("mov $0x0 , %ecx\n\t");
    __asm__("cpuid\n\t");
    __asm__("mov %%eax, %0\n\t":"=r" (ecx_max)); //gives maximum ECX value
    __asm__("mov %%ebx, %0\n\t":"=r" (ebx)); //extended feature flags
    __asm__("mov %%ecx, %0\n\t":"=r" (ecx)); //extended feature flags
    __asm__("mov %%edx, %0\n\t":"=r" (edx)); //extended feature flags
      
    // Intel-defined CPU features, CPUID level 0x00000007:0 (EBX), word 9
    {
      std::string features[] = { "fsgsbase",   /* RDFSBASE, WRFSBASE, RDGSBASE, WRGSBASE instructions*/
                                 "tsc_adjust", /* TSC adjustment MSR 0x3B */
                                 "sgx",        /* Software Guard Extensions */
                                 "bmi1",       /* 1st group bit manipulation extensions */
                                 "hle",        /* Hardware Lock Elision */
                                 "avx2",       /* AVX2 instructions */
                                 "",           /* "" FPU data pointer updated only on x87 exceptions */
                                 "smep",       /* Supervisor Mode Execution Protection */
                                 "bmi2",       /* 2nd group bit manipulation extensions */
                                 "erms",       /* Enhanced REP MOVSB/STOSB instructions */
                                 "invpcid",    /* Invalidate Processor Context ID */
                                 "rtm",        /* Restricted Transactional Memory */
                                 "cqm",        /* Cache QoS Monitoring */
                                 "",           /* "" Zero out FPU CS and FPU DS */
                                 "mpx",        /* Memory Protection Extension */
                                 "rdt_a",      /* Resource Director Technology Allocation */
                                 "avx512f",    /* AVX-512 Foundation */
                                 "avx512dq",   /* AVX-512 DQ (Double/Quad granular) Instructions */
                                 "rdseed",     /* RDSEED instruction */
                                 "adx",        /* ADCX and ADOX instructions */
                                 "smap",       /* Supervisor Mode Access Prevention */
                                 "avx512ifma", /* AVX-512 Integer Fused Multiply-Add instructions */
                                 "pcommit",
                                 "clflushopt", /* CLFLUSHOPT instruction */
                                 "clwb",       /* CLWB instruction */
                                 "intel_pt",   /* Intel Processor Trace */
                                 "avx512pf",   /* AVX-512 Prefetch */
                                 "avx512er",   /* AVX-512 Exponential and Reciprocal */
                                 "avx512cd",   /* AVX-512 Conflict Detection */
                                 "sha_ni",     /* SHA1/SHA256 Instruction Extensions */
                                 "avx512bw",   /* AVX-512 BW (Byte/Word granular) Instructions */
                                 "avx512vl"    /* AVX-512 VL (128/256 Vector Length) Extensions */
      };
      printf("INTEL 0x00000007:0 (EBX)");
      print_features(ebx, features, 32);
    }
      
    // Intel-defined CPU features, CPUID level 0x00000007:0 (ECX), word 16
    {
      std::string features[] = { "prefetchwt1", 
                                 "avx512vbmi",      /* AVX512 Vector Bit Manipulation instructions*/
                                 "umip",            /* User Mode Instruction Protection */
                                 "pku",             /* Protection Keys for Userspace */
                                 "ospke",           /* OS Protection Keys Enable */
                                 "waitpkg",         /* UMONITOR/UMWAIT/TPAUSE Instructions */
                                 "avx512vbmi2",     /* Additional AVX512 Vector Bit Manipulation Instructions */
                                 "cetss",
                                 "gfni",            /* Galois Field New Instructions */
                                 "vaes",            /* Vector AES */
                                 "vpclmulqdq",      /* Carry-Less Multiplication Double Quadword */
                                 "avx512vnni",      /* Vector Neural Network Instructions */
                                 "avx512bitalg",    /* Support for VPOPCNT[B,W] and VPSHUF-BITQMB instructions */
                                 "tme",             /* Intel Total Memory Encryption */
                                 "avx512vpopcntdq", /* POPCNT for vectors of DW/QW */
                                 "",                /* Reserved */
                                 "la57",            /* 5-level page tables */
                                 "",                /* Reserved */
                                 "",                /* Reserved */
                                 "",                /* Reserved */
                                 "",                /* Reserved */
                                 "",                /* Reserved */
                                 "rdpid",           /* RDPID instruction */
                                 "keylocker", 
                                 "bus_lock_detect", /* Bus Lock detect */
                                 "cldemote",        /* CLDEMOTE instruction */
                                 "",                /* Reserved */
                                 "movdiri",         /* MOVDIRI instruction */
                                 "movdir64b",       /* MOVDIR64B instruction */
                                 "enqcmd",          /* ENQCMD and ENQCMDS instructions */
                                 "sgx_lc",          /* Software Guard Extensions Launch Control */
                                 "pks"
      };
      printf("INTEL 0x00000007:0 (ECX)");
      print_features(ecx, features, 32);
    }
      
    // Intel-defined CPU features, CPUID level 0x00000007:0 (EDX), word 18
    {
      std::string features[] = { "",                   /* Reserved */
                                 "",                   /* Reserved */
                                 "avx5124vnniw",       /* AVX-512 Neural Network Instructions */
                                 "avx5124fmaps",       /* AVX-512 Multiply Accumulation Single precision */
                                 "fsrm",               /* Fast Short Rep Mov */
                                 "",                   /* Reserved */
                                 "",                   /* Reserved */
                                 "",                   /* Reserved */
                                 "avx512vp2intersect", /* AVX-512 Intersect for D/Q */
                                 "srbds",              /* "" SRBDS mitigation MSR available */
                                 "md_clear",           /* VERW clears CPU buffers */
                                 "",                   /* "" RTM transaction always aborts */
                                 "",                   /* Reserved */
                                 "",                   /* "" TSX_FORCE_ABORT */
                                 "serialize",          /* SERIALIZE instruction */
                                 "",                   /* "" This part has CPUs of more than one type */
                                 "tsxldtrk",           /* TSX Suspend Load Address Tracking */
                                 "",                   /* Reserved */
                                 "pconfig",            /* Intel PCONFIG */
                                 "arch_lbr",           /* Intel ARCH LBR */
                                 "cet_ibt",
                                 "",                   /* Reserved */
                                 "amx-bf16",           /* AMX BFLOAT16 Support */
                                 "avx512fp16",         /* AVX512 FP16 */
                                 "amx-tile",           /* AMX tile Support */
                                 "amx-int8",           /* AMX int8 Support */
                                 "ibrs ibpb",          /* "" Speculation Control (IBRS + IBPB) */
                                 "stibp",              /* "" Single Thread Indirect Branch Predictors */
                                 "flush_l1d",          /* Flush L1D cache */
                                 "arch_capabilities",  /* IA32_ARCH_CAPABILITIES MSR (Intel) */
                                 "",                   /* "" IA32_CORE_CAPABILITIES MSR */
                                 "ssbd"                /* "" Speculative Store Bypass Disable */
      };
      printf("INTEL 0x00000007:0 (EDX)");
      print_features(edx, features, 32);
    }

    printf("\necx_max:%d\n", ecx_max);
    
    if (ecx_max >= 1) {
      // EAX=7, ECX=1: Extended Features
      __asm__("mov $0x7 , %eax\n\t");
      __asm__("mov $0x1 , %ecx\n\t");
      __asm__("cpuid\n\t");
      __asm__("mov %%eax, %0\n\t":"=r" (eax)); //extended feature flags
      
      // Intel-defined CPU features, CPUID level 0x00000007:1 (EAX), word 12
      {
        std::string features[] = { "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "avx_vnni",   /* AVX VNNI instructions */
                                   "avx512bf16", /* AVX512 BFLOAT16 instructions */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   "",           /* Reserved */
                                   ""            /* Reserved */
        };
        printf("INTEL 0x00000007:1 (EAX)");
        print_features(eax, features, 32);
      }      
    } // ECX=1
  } // EAX=7

  // EAX=0: largest value that EAX can be set to before calling CPUID
  __asm__("mov $0x80000000, %eax\n\t");
  __asm__("cpuid\n\t");
  __asm__("mov %%eax, %0\n\t":"=r" (eax_max));

  printf("\neax_max:%d\n", eax_max);
  
  if (eax_max >= 1) {
    
    // EAX=80000001h: Processor Info and Feature Bits
    __asm__("mov $0x80000001 , %eax\n\t");
    __asm__("cpuid\n\t");
    __asm__("mov %%ecx, %0\n\t":"=r" (ecx)); //feature flags
    __asm__("mov %%edx, %0\n\t":"=r" (edx)); //feature flags
    
    // AMD-defined CPU features, CPUID level 0x80000001 (EDX), word 1
    // Don't duplicate feature flags which are redundant with Intel!
    {
      std::string features[] = { "",         /* Onboard FPU */
                                 "",         /* Virtual Mode Extensions */
                                 "",         /* Debugging Extensions */
                                 "",         /* Page Size Extensions */
                                 "",         /* Time Stamp Counter */
                                 "",         /* Model-Specific Registers */
                                 "",         /* Physical Address Extensions */
                                 "",         /* Machine Check Exception */
                                 "",         /* CMPXCHG8 instruction */
                                 "",         /* Onboard APIC */
                                 "",         /* Reserved */
                                 "syscall",  /* SYSCALL/SYSRET */
                                 "",         /* Memory Type Range Registers */
                                 "",         /* Page Global Enable */
                                 "",         /* Machine Check Architecture */
                                 "",         /* CMOV instructions (plus FCMOVcc, FCOMI with FPU) */
                                 "",         /* Page Attribute Table */
                                 "",         /* 36-bit PSEs */
                                 "",         /* Reserved */
                                 "mp",       /* MP Capable */
                                 "nx",       /* Execute Disable */
                                 "",         /* Reserved */
                                 "mmxext",   /* AMD MMX extensions */
                                 "",         /* Multimedia Extensions */
                                 "",         /* FXSAVE/FXRSTOR, CR4.OSFXSR */
                                 "fxsr_opt", /* FXSAVE/FXRSTOR optimizations */
                                 "pdpe1gb",  /* "pdpe1gb" GB pages */
                                 "rdtscp",   /* RDTSCP */
                                 "",         /* Reserved */
                                 "lm",       /* Long Mode (x86-64, 64-bit support) */
                                 "3dnowext", /* AMD 3DNow extensions */
                                 "3dnow"     /* 3DNow */
      };
      print_features(edx, features, 32);
    }
    
    // AMD-defined CPU features, CPUID level 0x80000001 (ECX), word 6
    {
      std::string features[] = { "lahf_lm",       /* LAHF/SAHF in long mode */
                                 "cmp_legacy",    /* If yes HyperThreading not valid */
                                 "svm",           /* Secure Virtual Machine */
                                 "extapic",       /* Extended APIC space */
                                 "cr8_legacy",    /* CR8 in 32-bit mode */
                                 "abm",           /* Advanced bit manipulation */
                                 "sse4a",         /* SSE-4A */
                                 "misalignsse",   /* Misaligned SSE mode */
                                 "3dnowprefetch", /* 3DNow prefetch instructions */
                                 "osvw",          /* OS Visible Workaround */
                                 "ibs",           /* Instruction Based Sampling */
                                 "xop",           /* extended AVX instructions */
                                 "skinit",        /* SKINIT/STGI instructions */
                                 "wdt",           /* Watchdog timer */
                                 "",              /* Reserved */
                                 "lwp",           /* Light Weight Profiling */
                                 "fma4",          /* 4 operands MAC instructions */
                                 "tce",           /* Translation Cache Extension */
                                 "",              /* Reserved */
                                 "nodeid_msr",    /* NodeId MSR */
                                 "",              /* Reserved */
                                 "tbm",           /* Trailing Bit Manipulations */
                                 "topoext",       /* Topology extensions CPUID leafs */
                                 "perfctr_core",  /* Core performance counter extensions */
                                 "perfctr_nb",    /* NB performance counter extensions */
                                 "",              /* Reserved */
                                 "bpext",         /* Data breakpoint extension */
                                 "ptsc",          /* Performance time-stamp counter */
                                 "perfctr_llc",   /* Last Level Cache performance counter extensions */
                                 "mwaitx",        /* MWAIT extension (MONITORX/MWAITX instructions) */
                                 "",              /* Reserved */
                                 ""               /* Reserved */
                                
      };
      print_features(ecx, features, 32);
    }
  } // EAX=1

  if (eax_max >=7) {      
    // EAX=7, ECX=0: Extended Features
    __asm__("mov $0x7 , %eax\n\t");
    __asm__("mov $0x0 , %ecx\n\t");
    __asm__("cpuid\n\t");
    __asm__("mov %%eax, %0\n\t":"=r" (ecx_max)); //gives maximum ECX value
    __asm__("mov %%ebx, %0\n\t":"=r" (ebx)); //extended feature flags
    __asm__("mov %%ecx, %0\n\t":"=r" (ecx)); //extended feature flags
    __asm__("mov %%edx, %0\n\t":"=r" (edx)); //extended feature flags
    
    //THIS IS NOT YET WORKING
    
    {
      std::string features[] = { "cxmmx",          /* Cyrix MMX extensions */
                                 "k6_mtrr",        /* AMD K6 nonstandard MTRRs */
                                 "cyrix_arr",      /* Cyrix ARRs (= MTRRs) */
                                 "centaur_mcr",    /* Centaur MCRs (= MTRRs) */
                                 "k8",             /* "" Opteron, Athlon64 */
                                 "",               /* "" Athlon */
                                 "",               /* "" P3 */
                                 "",               /* "" P4 */
                                 "constant_tsc",   /* TSC ticks at a constant rate */
                                 "up",             /* SMP kernel running on UP */
                                 "art",            /* Always running timer (ART) */
                                 "arch_perfmon",   /* Intel Architectural PerfMon */
                                 "pebs",           /* Precise-Event Based Sampling */
                                 "bts",            /* Branch Trace Store */
                                 "",               /* "" syscall in IA32 userspace */
                                 "",               /* "" sysenter in IA32 userspace */
                                 "rep_good",       /* REP microcode works well */
                                 "",               /* Reserved */
                                 "",               /* "" LFENCE synchronizes RDTSC */
                                 "acc_power",      /* AMD Accumulated Power Mechanism */
                                 "nopl",           /* The NOPL (0F 1F) instructions */
                                 "",               /* "" Always-present feature */
                                 "xtopology",      /* CPU topology enum extensions */
                                 "tsc_reliable",   /* TSC is known to be reliable */
                                 "nonstop_tsc",    /* TSC does not stop in C states */
                                 "cpuid",          /* CPU has CPUID instruction itself */
                                 "extd_apicid",    /* Extended APICID (8 bits) */
                                 "amd_dcm",        /* AMD multi-node processor */
                                 "aperfmperf",     /* P-State hardware coordination feedback capability (APERF/MPERF MSRs) */
                                 "rapl",           /* AMD/Hygon RAPL interface */
                                 "nonstop_tsc_s3", /* TSC doesn't stop in S3 state */
                                 "tsc_known_freq"  /* TSC has known frequency */                               
      };
      //print_features(eax, features, 32);
    }
  } // EAX=7
  
  printf("\n");
}

int main(){
  getVendorID();
  getProcInfo();
  getFeatures();    
  return 0;
}
