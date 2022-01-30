#include <cstdio>
#include <cstring>
#include <string>

int main(){
  int a[4];
  for(int i=0; i<4; ++i)
    a[i] = 0;
  
  // EAX=0: Highest Function Parameter and Manufacturer ID
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

  // EAX=1: Processor Info and Feature Bits
  __asm__("mov $0x1 , %eax\n\t");
  __asm__("cpuid\n\t");
  __asm__("mov %%eax, %0\n\t":"=r" (a[0])); //gives model and family
  __asm__("mov %%ebx, %0\n\t":"=r" (a[1])); //gives additional feature info
  __asm__("mov %%ecx, %0\n\t":"=r" (a[2])); //feature flags
  __asm__("mov %%edx, %0\n\t":"=r" (a[3])); //feature flags

  int stepping = a[0]>>0 & 0xF;
  int model    = a[0]>>4 & 0xF;
  int family   = a[0]>>8 & 0xF;
  if(family == 6 || family == 15)
    model += (a[0]>>16 & 0xF)<<4;

  printf ("cpu family      : %d\n", family);
  printf ("model           : %d\n", model);
  printf ("stepping        : %d\n", stepping);

  // CPU flags
  printf ("flags           : ");

  // Features in EDX register
  std::string edx_feature[] = {
    "fpu", "vme", "de", "pse",
    "tsc", "msr", "pae", "mce",
    "cx8", "apic" , "", "sep",
    "mtrr", "pge", "mca", "cmov",
    "pat", "pse36", "psn", "clflush",
    "", "dts", "acpi", "mmx",
    "fxsr", "sse", "sse2", "ss",
    "htt", "tm", "ia64", "pbe" };

  for (int i=0; i<32; ++i)
    printf ("%s", (a[3]>>i & 0x1) && !edx_feature[i].empty() ? (edx_feature[i]+" ").c_str() : "");
  
  // Features in ECX register
  std::string ecx_feature[] = {
    "sse3", "pclmulqdq", "dtes64", "monitor",
    "ds-cpl", "vmx", "smx", "est",
    "tm2", "ssse3", "cnxt-id", "sdbg",
    "fma", "cx16", "xtpr", "pdcm",
    "", "pcid", "dca", "sse4_1",
    "sse4_2", "x2apic", "movbe", "popcnt",
    "tsc-deadline", "aes", "xsave", "osxsave",
    "avx", "f16c", "rdrnd", "hypervisor"
  };

  for (int i=0; i<32; ++i)
    printf ("%s", (a[2]>>i & 0x1) && !ecx_feature[i].empty() ? (ecx_feature[i]+" ").c_str() : "");

  // EAX=7: Extended Features
  __asm__("mov $0x7 , %eax\n\t");
  __asm__("mov $0x0 , %ecx\n\t");
  __asm__("cpuid\n\t");
  __asm__("mov %%eax, %0\n\t":"=r" (a[0])); //gives maximum ECX value
  __asm__("mov %%ebx, %0\n\t":"=r" (a[1])); //extended feature flags
  __asm__("mov %%ecx, %0\n\t":"=r" (a[2])); //extended feature flags
  __asm__("mov %%edx, %0\n\t":"=r" (a[3])); //extended feature flags
  
  // Extended features in EBX register
  std::string ebx_extended_feature[] = {
    "fsgsbase", "", "sgx", "bmi1",
    "hle", "avx2", "", "smep",
    "bmi2", "erms", "invpcid", "rtm",
    "pqm", "", "mpx", "pqe",
    "avx512f", "avx512dq", "rdseed", "adx",
    "smap", "avx512ifma", "pcommit", "clflushopt",
    "clwb", "intelpt", "avx512pf", "avx512er",
    "avx512cd", "sha", "avx512bw", "avx512vl"
  };
  
  for (int i=0; i<32; ++i)
    printf ("%s", (a[1]>>i & 0x1) && !ebx_extended_feature[i].empty() ? (ebx_extended_feature[i]+" ").c_str() : "");

  // Extended features in ECX register
  std::string ecx_extended_feature[] = {
    "prefetchwt1", "avx512vbmi", "umip", "pku",
    "ospke", "waitpkg", "avx512vbmi2", "cetss",
    "gfni", "vaes", "vpclmulqdq", "avx512vnni",
    "avx512bitalg", "TMEEN", "avx512vpopcntdq", "",
    "", "", "", "",
    "", "", "rdpid", "keylocker",
    "", "cldemote", "", "movdiri",
    "movdir64b", "enqcmd", "sgx_lc", "pks"
  };

  for (int i=0; i<32; ++i)
    printf ("%s", (a[2]>>i & 0x1) && !ecx_extended_feature[i].empty() ? (ecx_extended_feature[i]+" ").c_str() : "");

  // Extended features in EDX register
  std::string edx_extended_feature[] = {
    "", "", "avx5124vnniw", "avx5124fmaps",
    "fsrm", "", "", "",
    "avx512vp2intersect", "SRBDS_CTRL", "md_clear", "",
    "", "tsx_force_abort", "serialize", "hybrid",
    "tsxldtrk", "", "pconfig", "lbr",
    "cet_ibt", "", "amx-bf16", "avx512fp16",
    "amx-tile", "amx-int8", "IBRS_IBPB", "stibp",
    "L1D_FLUSH", "IA32_ARCH_CAPABILITIES", "IA32_CORE_CAPABILITIES", "ssbd"
  };

  for (int i=0; i<32; ++i)
    printf ("%s", (a[3]>>i & 0x1) && !edx_extended_feature[i].empty() ? (edx_extended_feature[i]+" ").c_str() : "");
  
  printf("\n");
  return  0;
}
