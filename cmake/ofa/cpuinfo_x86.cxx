#include <cstdio>
#include <cstring>
#include <string>

#define print_features(reg,features,n)                                  \
  for (int i=0; i<n; ++i)                                               \
    printf ("%s", (reg>>i & 0x1) && !features[i].empty()                \
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
  int32_t a[3], eax,ebx,ecx,edx;
  for(int i=0; i<3; ++i)
    a[i] = 0;

  // CPU flags
  printf ("flags           : ");
  
  // EAX=0: Vendor ID
  __asm__("mov $0x0, %eax\n\t");
  __asm__("cpuid\n\t");
  __asm__("mov %%eax, %0\n\t":"=r" (eax));
  __asm__("mov %%ebx, %0\n\t":"=r" (a[0]));
  __asm__("mov %%edx, %0\n\t":"=r" (a[1]));
  __asm__("mov %%ecx, %0\n\t":"=r" (a[2]));

  char vendorID[13]; vendorID[12] = 0;
  memcpy(&vendorID[0],&a[0],4);
  memcpy(&vendorID[4],&a[1],4);
  memcpy(&vendorID[8],&a[2],4);
  
  if (strcmp(vendorID, "GenuineIntel") == 0) {

    if (eax >= 1) {
      
      // EAX=1: Processor Info and Feature Bits
      __asm__("mov $0x1 , %eax\n\t");
      __asm__("cpuid\n\t");
      __asm__("mov %%ecx, %0\n\t":"=r" (ecx)); //feature flags
      __asm__("mov %%edx, %0\n\t":"=r" (edx)); //feature flags   
      
      // Features in EDX register
      {
        std::string features[] = {
          "fpu", "vme", "de", "pse",
          "tsc", "msr", "pae", "mce",
          "cx8", "apic" , "", "sep",
          "mtrr", "pge", "mca", "cmov",
          "pat", "pse36", "psn", "clflush",
          "", "dts", "acpi", "mmx",
          "fxsr", "sse", "sse2", "ss",
          "htt", "tm", "ia64", "pbe" };
        print_features(edx, features, 32);
      }
      
      // Features in ECX register
      {
        std::string features[] = {
          "sse3", "pclmulqdq", "dtes64", "monitor",
          "ds-cpl", "vmx", "smx", "est",
          "tm2", "ssse3", "cnxt-id", "sdbg",
          "fma", "cx16", "xtpr", "pdcm",
          "", "pcid", "dca", "sse4_1",
          "sse4_2", "x2apic", "movbe", "popcnt",
          "tsc-deadline", "aes", "xsave", "osxsave",
          "avx", "f16c", "rdrnd", "hypervisor"
        };
        print_features(ecx, features, 32);
      }
    }
    
    if (eax >=7) {      
      // EAX=7, ECX=0: Extended Features
      __asm__("mov $0x7 , %eax\n\t");
      __asm__("mov $0x0 , %ecx\n\t");
      __asm__("cpuid\n\t");
      __asm__("mov %%eax, %0\n\t":"=r" (eax)); //gives maximum ECX value
      __asm__("mov %%ebx, %0\n\t":"=r" (ebx)); //extended feature flags
      __asm__("mov %%ecx, %0\n\t":"=r" (ecx)); //extended feature flags
      __asm__("mov %%edx, %0\n\t":"=r" (edx)); //extended feature flags
      
      // Extended features in EBX register
      {
        std::string features[] = {
          "fsgsbase", "", "sgx", "bmi1",
          "hle", "avx2", "", "smep",
          "bmi2", "erms", "invpcid", "rtm",
          "pqm", "", "mpx", "pqe",
          "avx512f", "avx512dq", "rdseed", "adx",
          "smap", "avx512ifma", "pcommit", "clflushopt",
          "clwb", "intelpt", "avx512pf", "avx512er",
          "avx512cd", "sha", "avx512bw", "avx512vl"
        };
        print_features(ebx, features, 32);
      }
      
      // Extended features in ECX register
      {
        std::string features[] = {
          "prefetchwt1", "avx512vbmi", "umip", "pku",
          "ospke", "waitpkg", "avx512vbmi2", "cetss",
          "gfni", "vaes", "vpclmulqdq", "avx512vnni",
          "avx512bitalg", "TMEEN", "avx512vpopcntdq", "",
          "", "", "", "",
          "", "", "rdpid", "keylocker",
          "", "cldemote", "", "movdiri",
          "movdir64b", "enqcmd", "sgx_lc", "pks"
        };
        print_features(ecx, features, 32);
      }
      
      // Extended features in EDX register
      {
        std::string features[] = {
          "", "", "avx5124vnniw", "avx5124fmaps",
          "fsrm", "", "", "",
          "avx512vp2intersect", "SRBDS_CTRL", "md_clear", "",
          "", "tsx_force_abort", "serialize", "hybrid",
          "tsxldtrk", "", "pconfig", "lbr",
          "cet_ibt", "", "amx-bf16", "avx512fp16",
          "amx-tile", "amx-int8", "IBRS_IBPB", "stibp",
          "L1D_FLUSH", "IA32_ARCH_CAPABILITIES", "IA32_CORE_CAPABILITIES", "ssbd"
        };
        print_features(edx, features, 32);
      }

      if (eax >= 1) {

        // Extended features in EAX register
        {
          std::string features[] = {
            "", "", "", "",
            "", "avx512bf16", "", "",
            "", "", "", "",
            "", "", "", "",
            "", "", "", "",
            "", "", "", "",
            "", "", "", "",
            "", "", "", ""
          };
          print_features(eax, features, 32);
        }
      }
    }
  }
  else if (strcmp(vendorID, "AuthenticAMD") == 0) {

    // // EAX=7, ECX=1: Extended Features
    // __asm__("mov $0x7 , %eax\n\t");
    // __asm__("mov $0x1 , %ecx\n\t");
    // __asm__("cpuid\n\t");
    // __asm__("mov %%eax, %0\n\t":"=r" (eax)); //extended feature flags
    
  }
  
  printf("\n");
}

int main(){
  getVendorID();
  getProcInfo();
  getFeatures();    
  return 0;
}
