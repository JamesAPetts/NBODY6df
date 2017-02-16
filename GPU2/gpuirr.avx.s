	.file	"gpuirr.avx.cpp"
	.text
	.p2align 4,,15
	.type	_ZL11gpuirr_openii._omp_fn.1, @function
_ZL11gpuirr_openii._omp_fn.1:
.LFB1269:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	call	omp_get_num_threads
	movl	%eax, _ZL11num_threads(%rip)
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1269:
	.size	_ZL11gpuirr_openii._omp_fn.1, .-_ZL11gpuirr_openii._omp_fn.1
	.p2align 4,,15
	.type	_ZL11gpuirr_firriPdS_, @function
_ZL11gpuirr_firriPdS_:
.LFB1257:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movslq	%edi, %rdi
	leaq	(%rdi,%rdi,2), %rax
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	imulq	$10016, %rdi, %r8
	salq	$5, %rax
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	andq	$-32, %rsp
	subq	$456, %rsp
	movq	%rdx, 432(%rsp)
	movq	_ZL4ptcl(%rip), %rdx
	vmovapd	_ZL8vec_tnow(%rip), %xmm0
	addq	_ZL4list(%rip), %r8
	vmovapd	%xmm0, 440(%rsp)
	addq	%rdx, %rax
	vmovddup	80(%rax), %xmm4
	vsubpd	%xmm4, %xmm0, %xmm4
	movl	12(%r8), %ecx
	vcvtpd2psx	%xmm4, %xmm4
	vshufps	$0, %xmm4, %xmm4, %xmm0
	vmovaps	%xmm0, %xmm0
	vshufps	$85, %xmm4, %xmm4, %xmm4
	vinsertf128	$0x1, %xmm4, %ymm0, %ymm4
	vmovaps	16(%rax), %xmm0
	vaddps	%ymm4, %ymm4, %ymm7
	testl	%ecx, %ecx
	vmovaps	%xmm0, %xmm2
	vinsertf128	$0x1, %xmm0, %ymm2, %ymm2
	vmovaps	32(%rax), %xmm0
	vmovaps	%xmm0, %xmm1
	vinsertf128	$0x1, %xmm0, %ymm1, %ymm0
	vmovaps	48(%rax), %xmm1
	vmovaps	%xmm1, %xmm5
	vinsertf128	$0x1, %xmm1, %ymm5, %ymm5
	vmovaps	64(%rax), %xmm1
	vmovaps	%xmm1, %xmm6
	vinsertf128	$0x1, %xmm1, %ymm6, %ymm6
	vmovaps	(%rax), %xmm1
	vmovaps	%xmm1, %xmm3
	vinsertf128	$0x1, %xmm1, %ymm3, %ymm3
	vmulps	%ymm6, %ymm4, %ymm1
	vaddps	%ymm1, %ymm5, %ymm1
	vmulps	%ymm1, %ymm4, %ymm1
	vaddps	%ymm1, %ymm0, %ymm1
	vmulps	%ymm1, %ymm4, %ymm1
	vaddps	%ymm1, %ymm2, %ymm2
	vmovaps	.LC0(%rip), %ymm1
	vmulps	%ymm1, %ymm4, %ymm4
	vmulps	%ymm6, %ymm4, %ymm6
	vaddps	%ymm6, %ymm5, %ymm5
	vshufps	$0, %ymm2, %ymm2, %ymm6
	vmovaps	%ymm6, 72(%rsp)
	vmulps	%ymm5, %ymm7, %ymm4
	vaddps	%ymm4, %ymm0, %ymm0
	vshufps	$85, %ymm3, %ymm3, %ymm5
	vshufps	$0, %ymm3, %ymm3, %ymm4
	vshufps	$85, %ymm2, %ymm2, %ymm7
	vshufps	$170, %ymm3, %ymm3, %ymm3
	vshufps	$170, %ymm2, %ymm2, %ymm2
	vmovaps	%ymm3, 104(%rsp)
	vmovaps	%ymm2, 8(%rsp)
	vshufps	$85, %ymm0, %ymm0, %ymm3
	vshufps	$0, %ymm0, %ymm0, %ymm2
	vmovaps	%ymm4, 168(%rsp)
	vshufps	$170, %ymm0, %ymm0, %ymm0
	vmovaps	%ymm5, 136(%rsp)
	vmovaps	%ymm7, 40(%rsp)
	vmovaps	%ymm2, -24(%rsp)
	vmovaps	%ymm3, -56(%rsp)
	vmovaps	%ymm0, -88(%rsp)
	jle	.L6
	.cfi_offset 3, -56
	.cfi_offset 12, -48
	.cfi_offset 13, -40
	.cfi_offset 14, -32
	.cfi_offset 15, -24
	leal	-1(%rcx), %r9d
	movl	%r9d, %eax
	leaq	4(,%rax,4), %rdi
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L4:
	movslq	16(%r8,%rax), %rcx
	addq	$4, %rax
	leaq	(%rcx,%rcx,2), %rcx
	salq	$5, %rcx
	addq	%rdx, %rcx
	cmpq	%rdi, %rax
	prefetcht0	(%rcx)
	prefetcht0	64(%rcx)
	jne	.L4
	vxorps	%xmm4, %xmm4, %xmm4
	leaq	4(%r8), %r15
	shrl	$3, %r9d
	vmovaps	%ymm4, 360(%rsp)
	addq	$1, %r9
	xorl	%eax, %eax
	vmovaps	%ymm4, 296(%rsp)
	salq	$5, %r9
	vmovaps	%ymm4, 328(%rsp)
	leaq	12(%r8), %r14
	vmovaps	%ymm4, 200(%rsp)
	leaq	8(%r8), %r13
	vmovaps	%ymm4, 232(%rsp)
	leaq	20(%r8), %r12
	vmovaps	%ymm4, 264(%rsp)
	leaq	16(%r8), %rbx
	leaq	28(%r8), %r11
	leaq	24(%r8), %r10
	.p2align 4,,10
	.p2align 3
.L5:
	movslq	16(%r8,%rax), %rdi
	movslq	16(%r15,%rax), %rcx
	vmovapd	440(%rsp), %xmm5
	leaq	(%rdi,%rdi,2), %rdi
	leaq	(%rcx,%rcx,2), %rcx
	salq	$5, %rdi
	addq	%rdx, %rdi
	salq	$5, %rcx
	vmovapd	80(%rdi), %xmm0
	addq	%rdx, %rcx
	vmovaps	64(%rdi), %xmm6
	vmovhpd	80(%rcx), %xmm0, %xmm0
	vmovaps	16(%rdi), %xmm8
	vsubpd	%xmm0, %xmm5, %xmm0
	vmovaps	48(%rdi), %xmm5
	vmovaps	%xmm6, %xmm6
	vinsertf128	$0x1, 64(%rcx), %ymm6, %ymm6
	vmovaps	(%rdi), %xmm12
	vmovaps	%xmm8, %xmm8
	vmovaps	%xmm5, %xmm5
	vinsertf128	$0x1, 16(%rcx), %ymm8, %ymm8
	vcvtpd2psx	%xmm0, %xmm0
	vshufps	$0, %xmm0, %xmm0, %xmm2
	vmovaps	%xmm2, %xmm2
	vinsertf128	$0x1, 48(%rcx), %ymm5, %ymm5
	vshufps	$85, %xmm0, %xmm0, %xmm0
	vmovaps	%xmm12, %xmm12
	vinsertf128	$0x1, %xmm0, %ymm2, %ymm0
	vmovaps	32(%rdi), %xmm2
	vmulps	%ymm6, %ymm0, %ymm3
	vaddps	%ymm3, %ymm5, %ymm3
	movslq	16(%r13,%rax), %rdi
	vinsertf128	$0x1, (%rcx), %ymm12, %ymm12
	vmovaps	%xmm2, %xmm2
	vinsertf128	$0x1, 32(%rcx), %ymm2, %ymm2
	movslq	16(%r14,%rax), %rcx
	vmulps	%ymm3, %ymm0, %ymm3
	vaddps	%ymm3, %ymm2, %ymm3
	leaq	(%rdi,%rdi,2), %rdi
	salq	$5, %rdi
	vmulps	%ymm3, %ymm0, %ymm3
	vaddps	%ymm3, %ymm8, %ymm8
	leaq	(%rcx,%rcx,2), %rcx
	vaddps	%ymm0, %ymm0, %ymm3
	vmulps	%ymm1, %ymm0, %ymm0
	vmulps	%ymm6, %ymm0, %ymm0
	vaddps	%ymm0, %ymm5, %ymm0
	addq	%rdx, %rdi
	salq	$5, %rcx
	addq	%rdx, %rcx
	vmovaps	48(%rdi), %xmm7
	vmovapd	440(%rsp), %xmm6
	vmulps	%ymm0, %ymm3, %ymm0
	vaddps	%ymm0, %ymm2, %ymm4
	vmovaps	%xmm7, %xmm7
	vmovapd	80(%rdi), %xmm0
	vinsertf128	$0x1, 48(%rcx), %ymm7, %ymm7
	vmovaps	64(%rdi), %xmm11
	vmovhpd	80(%rcx), %xmm0, %xmm0
	vmovaps	32(%rdi), %xmm3
	vsubpd	%xmm0, %xmm6, %xmm0
	vmovaps	%xmm11, %xmm11
	vmovaps	16(%rdi), %xmm6
	vinsertf128	$0x1, 64(%rcx), %ymm11, %ymm11
	vmovaps	(%rdi), %xmm10
	vmovaps	%xmm3, %xmm3
	movslq	16(%rbx,%rax), %rdi
	vinsertf128	$0x1, 32(%rcx), %ymm3, %ymm3
	vmovaps	%xmm6, %xmm6
	vcvtpd2psx	%xmm0, %xmm0
	vshufps	$0, %xmm0, %xmm0, %xmm2
	vmovaps	%xmm2, %xmm2
	vshufps	$85, %xmm0, %xmm0, %xmm0
	vinsertf128	$0x1, %xmm0, %ymm2, %ymm0
	vmulps	%ymm11, %ymm0, %ymm2
	vaddps	%ymm2, %ymm7, %ymm2
	vmovaps	%xmm10, %xmm10
	vinsertf128	$0x1, 16(%rcx), %ymm6, %ymm6
	vinsertf128	$0x1, (%rcx), %ymm10, %ymm10
	movslq	16(%r12,%rax), %rcx
	vaddps	%ymm0, %ymm0, %ymm5
	leaq	(%rdi,%rdi,2), %rdi
	vmulps	%ymm2, %ymm0, %ymm2
	vaddps	%ymm2, %ymm3, %ymm2
	salq	$5, %rdi
	addq	%rdx, %rdi
	leaq	(%rcx,%rcx,2), %rcx
	vmovaps	48(%rdi), %xmm9
	vmulps	%ymm2, %ymm0, %ymm2
	vmulps	%ymm1, %ymm0, %ymm0
	vmulps	%ymm11, %ymm0, %ymm0
	vaddps	%ymm0, %ymm7, %ymm0
	salq	$5, %rcx
	vmovaps	64(%rdi), %xmm13
	vaddps	%ymm2, %ymm6, %ymm6
	addq	%rdx, %rcx
	vmovaps	%xmm9, %xmm9
	vmovaps	%xmm13, %xmm13
	vinsertf128	$0x1, 48(%rcx), %ymm9, %ymm9
	vmovaps	(%rdi), %xmm11
	vmulps	%ymm0, %ymm5, %ymm0
	vaddps	%ymm0, %ymm3, %ymm2
	vmovapd	80(%rdi), %xmm0
	vinsertf128	$0x1, 64(%rcx), %ymm13, %ymm13
	vmovaps	32(%rdi), %xmm5
	vmovaps	%xmm11, %xmm11
	vmovapd	440(%rsp), %xmm7
	vmovhpd	80(%rcx), %xmm0, %xmm0
	vinsertf128	$0x1, (%rcx), %ymm11, %ymm11
	vmovaps	%xmm5, %xmm5
	vsubpd	%xmm0, %xmm7, %xmm0
	vinsertf128	$0x1, 32(%rcx), %ymm5, %ymm5
	vmovaps	16(%rdi), %xmm7
	vmovaps	%xmm7, %xmm7
	vcvtpd2psx	%xmm0, %xmm0
	vshufps	$0, %xmm0, %xmm0, %xmm3
	vmovaps	%xmm3, %xmm3
	vshufps	$85, %xmm0, %xmm0, %xmm0
	vinsertf128	$0x1, %xmm0, %ymm3, %ymm0
	vmulps	%ymm13, %ymm0, %ymm3
	vaddps	%ymm3, %ymm9, %ymm3
	vinsertf128	$0x1, 16(%rcx), %ymm7, %ymm7
	vaddps	%ymm0, %ymm0, %ymm14
	vmulps	%ymm3, %ymm0, %ymm3
	vaddps	%ymm3, %ymm5, %ymm3
	vmulps	%ymm3, %ymm0, %ymm3
	vmulps	%ymm1, %ymm0, %ymm0
	vmulps	%ymm13, %ymm0, %ymm0
	vaddps	%ymm0, %ymm9, %ymm0
	vaddps	%ymm3, %ymm7, %ymm7
	vmulps	%ymm0, %ymm14, %ymm0
	vaddps	%ymm0, %ymm5, %ymm0
	vmovaps	%ymm0, -120(%rsp)
	movslq	16(%r10,%rax), %rdi
	movslq	16(%r11,%rax), %rcx
	addq	$32, %rax
	vmovapd	440(%rsp), %xmm3
	leaq	(%rdi,%rdi,2), %rdi
	leaq	(%rcx,%rcx,2), %rcx
	salq	$5, %rdi
	addq	%rdx, %rdi
	salq	$5, %rcx
	vmovapd	80(%rdi), %xmm0
	addq	%rdx, %rcx
	cmpq	%r9, %rax
	vmovaps	48(%rdi), %xmm14
	vmovhpd	80(%rcx), %xmm0, %xmm0
	vmovaps	64(%rdi), %xmm15
	vsubpd	%xmm0, %xmm3, %xmm0
	vmovaps	%xmm14, %xmm14
	vmovaps	32(%rdi), %xmm13
	vmovaps	%xmm15, %xmm15
	vinsertf128	$0x1, 48(%rcx), %ymm14, %ymm14
	vinsertf128	$0x1, 64(%rcx), %ymm15, %ymm15
	vmovaps	%xmm13, %xmm13
	vcvtpd2psx	%xmm0, %xmm0
	vshufps	$85, %xmm0, %xmm0, %xmm9
	vshufps	$0, %xmm0, %xmm0, %xmm5
	vmovaps	%xmm5, %xmm0
	vmovaps	16(%rdi), %xmm5
	vinsertf128	$0x1, 32(%rcx), %ymm13, %ymm13
	vinsertf128	$0x1, %xmm9, %ymm0, %ymm0
	vmovaps	(%rdi), %xmm9
	vmovaps	%xmm5, %xmm5
	vinsertf128	$0x1, 16(%rcx), %ymm5, %ymm3
	vmulps	%ymm15, %ymm0, %ymm5
	vaddps	%ymm5, %ymm14, %ymm5
	vmovaps	%xmm9, %xmm9
	vinsertf128	$0x1, (%rcx), %ymm9, %ymm9
	vmulps	%ymm5, %ymm0, %ymm5
	vaddps	%ymm5, %ymm13, %ymm5
	vmulps	%ymm5, %ymm0, %ymm5
	vaddps	%ymm3, %ymm5, %ymm5
	vaddps	%ymm0, %ymm0, %ymm3
	vmulps	%ymm1, %ymm0, %ymm0
	vmulps	%ymm15, %ymm0, %ymm0
	vaddps	%ymm0, %ymm14, %ymm0
	vunpcklps	%ymm9, %ymm10, %ymm14
	vunpckhps	%ymm9, %ymm10, %ymm9
	vmulps	%ymm3, %ymm0, %ymm0
	vaddps	%ymm0, %ymm13, %ymm13
	vunpcklps	%ymm11, %ymm12, %ymm0
	vunpckhps	%ymm11, %ymm12, %ymm11
	vunpcklps	%ymm14, %ymm0, %ymm15
	vunpcklps	%ymm9, %ymm11, %ymm10
	vunpckhps	%ymm14, %ymm0, %ymm14
	vunpckhps	%ymm9, %ymm11, %ymm11
	vunpcklps	%ymm7, %ymm8, %ymm0
	vunpcklps	%ymm5, %ymm6, %ymm9
	vunpcklps	%ymm9, %ymm0, %ymm12
	vsubps	168(%rsp), %ymm15, %ymm15
	vunpckhps	%ymm9, %ymm0, %ymm9
	vsubps	72(%rsp), %ymm12, %ymm12
	vunpckhps	%ymm5, %ymm6, %ymm5
	vsubps	40(%rsp), %ymm9, %ymm9
	vunpckhps	%ymm7, %ymm8, %ymm7
	vunpcklps	%ymm5, %ymm7, %ymm7
	vsubps	104(%rsp), %ymm10, %ymm10
	vunpcklps	-120(%rsp), %ymm4, %ymm0
	vunpcklps	%ymm13, %ymm2, %ymm5
	vaddps	%ymm12, %ymm15, %ymm15
	vunpcklps	%ymm5, %ymm0, %ymm6
	vunpckhps	%ymm13, %ymm2, %ymm13
	vsubps	8(%rsp), %ymm7, %ymm7
	vunpckhps	%ymm5, %ymm0, %ymm5
	vunpckhps	-120(%rsp), %ymm4, %ymm3
	vsubps	136(%rsp), %ymm14, %ymm14
	vunpcklps	%ymm13, %ymm3, %ymm3
	vmulps	%ymm15, %ymm15, %ymm2
	vsubps	-24(%rsp), %ymm6, %ymm6
	vaddps	%ymm7, %ymm10, %ymm7
	vaddps	%ymm9, %ymm14, %ymm14
	vsubps	-56(%rsp), %ymm5, %ymm5
	vsubps	-88(%rsp), %ymm3, %ymm3
	vmulps	%ymm14, %ymm14, %ymm0
	vaddps	%ymm0, %ymm2, %ymm2
	vmulps	%ymm7, %ymm7, %ymm0
	vaddps	%ymm0, %ymm2, %ymm2
	vrsqrtps	%ymm2, %ymm0
	vmulps	%ymm0, %ymm2, %ymm2
	vmulps	%ymm2, %ymm0, %ymm2
	vaddps	.LC2(%rip), %ymm2, %ymm9
	vmulps	.LC1(%rip), %ymm0, %ymm4
	vmulps	%ymm5, %ymm14, %ymm2
	vmulps	%ymm9, %ymm4, %ymm9
	vmulps	%ymm6, %ymm15, %ymm4
	vaddps	%ymm2, %ymm4, %ymm4
	vmulps	%ymm3, %ymm7, %ymm2
	vmulps	%ymm9, %ymm9, %ymm0
	vmulps	%ymm9, %ymm11, %ymm11
	vaddps	%ymm2, %ymm4, %ymm4
	vmulps	.LC2(%rip), %ymm0, %ymm2
	vmulps	%ymm11, %ymm0, %ymm0
	vmulps	%ymm2, %ymm4, %ymm2
	vmulps	%ymm0, %ymm15, %ymm4
	vmulps	%ymm2, %ymm15, %ymm15
	vaddps	328(%rsp), %ymm4, %ymm4
	vaddps	%ymm15, %ymm6, %ymm6
	vmovaps	%ymm4, 328(%rsp)
	vmulps	%ymm0, %ymm14, %ymm4
	vaddps	296(%rsp), %ymm4, %ymm4
	vmulps	%ymm2, %ymm14, %ymm14
	vmulps	%ymm2, %ymm7, %ymm2
	vaddps	%ymm14, %ymm5, %ymm5
	vaddps	%ymm2, %ymm3, %ymm3
	vmulps	%ymm6, %ymm0, %ymm6
	vmovaps	%ymm4, 296(%rsp)
	vmulps	%ymm0, %ymm7, %ymm4
	vaddps	360(%rsp), %ymm4, %ymm4
	vmulps	%ymm5, %ymm0, %ymm5
	vmulps	%ymm3, %ymm0, %ymm3
	vmovaps	%ymm4, 360(%rsp)
	vaddps	264(%rsp), %ymm6, %ymm6
	vaddps	232(%rsp), %ymm5, %ymm5
	vaddps	200(%rsp), %ymm3, %ymm3
	vmovaps	%ymm6, 264(%rsp)
	vmovaps	%ymm5, 232(%rsp)
	vmovaps	%ymm3, 200(%rsp)
	jne	.L5
	vmovaps	328(%rsp), %ymm5
	vmovaps	296(%rsp), %ymm6
	vmovaps	360(%rsp), %ymm7
.L3:
	vcvtps2pd	328(%rsp), %ymm1
	vextractf128	$0x1, %ymm5, %xmm0
	vcvtps2pd	%xmm0, %ymm0
	vaddpd	%ymm0, %ymm1, %ymm0
	vcvtps2pd	296(%rsp), %ymm2
	vmovaps	264(%rsp), %ymm4
	movq	432(%rsp), %rax
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm1
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm1, %xmm1
	vextractf128	$0x1, %ymm6, %xmm0
	vcvtps2pd	%xmm0, %ymm0
	vmovaps	232(%rsp), %ymm6
	vaddpd	%ymm0, %ymm2, %ymm0
	vcvtps2pd	360(%rsp), %ymm2
	vmovlpd	%xmm1, (%rsi)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm3
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm3, %xmm3
	vextractf128	$0x1, %ymm7, %xmm0
	vcvtps2pd	%xmm0, %ymm0
	vmovaps	200(%rsp), %ymm7
	vaddpd	%ymm0, %ymm2, %ymm0
	vcvtps2pd	264(%rsp), %ymm2
	vmovlpd	%xmm3, 8(%rsi)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm5
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm5, %xmm5
	vextractf128	$0x1, %ymm4, %xmm0
	vcvtps2pd	%xmm0, %ymm0
	vaddpd	%ymm0, %ymm2, %ymm0
	vcvtps2pd	232(%rsp), %ymm2
	vmovlpd	%xmm5, 16(%rsi)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm4
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm4, %xmm4
	vextractf128	$0x1, %ymm6, %xmm0
	vcvtps2pd	%xmm0, %ymm0
	vcvtps2pd	200(%rsp), %ymm6
	vaddpd	%ymm0, %ymm2, %ymm0
	vmovlpd	%xmm4, (%rax)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm2
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm2, %xmm2
	vextractf128	$0x1, %ymm7, %xmm0
	vcvtps2pd	%xmm0, %ymm0
	vaddpd	%ymm0, %ymm6, %ymm0
	vmovlpd	%xmm2, 8(%rax)
	vhaddpd	%ymm0, %ymm0, %ymm0
	vmovapd	%xmm0, %xmm6
	vextractf128	$0x1, %ymm0, %xmm0
	vaddpd	%xmm0, %xmm6, %xmm0
	vmovlpd	%xmm0, 16(%rax)
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	vzeroupper
	ret
.L6:
	.cfi_restore_state
	vxorps	%xmm4, %xmm4, %xmm4
	vmovaps	%ymm4, %ymm5
	vmovaps	%ymm4, %ymm6
	vmovaps	%ymm4, 360(%rsp)
	vmovaps	%ymm4, %ymm7
	vmovaps	%ymm4, 296(%rsp)
	vmovaps	%ymm4, 328(%rsp)
	vmovaps	%ymm4, 200(%rsp)
	vmovaps	%ymm4, 232(%rsp)
	vmovaps	%ymm4, 264(%rsp)
	jmp	.L3
	.cfi_endproc
.LFE1257:
	.size	_ZL11gpuirr_firriPdS_, .-_ZL11gpuirr_firriPdS_
	.p2align 4,,15
	.type	_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.0, @function
_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.0:
.LFB1268:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movl	$1, %ecx
	movl	$1, %edx
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	movq	%rdi, %rbx
	.cfi_offset 3, -56
	.cfi_offset 12, -48
	.cfi_offset 13, -40
	.cfi_offset 14, -32
	.cfi_offset 15, -24
	andq	$-32, %rsp
	subq	$96, %rsp
	movslq	24(%rdi), %rsi
	leaq	64(%rsp), %r9
	xorl	%edi, %edi
	leaq	32(%rsp), %r8
	call	GOMP_loop_guided_start
	xorl	%r12d, %r12d
	testb	%al, %al
	je	.L11
	movq	_ZL4list(%rip), %r13
	.p2align 4,,10
	.p2align 3
.L13:
	movq	32(%rsp), %rax
	movq	(%rbx), %rcx
	movq	%rax, 24(%rsp)
	movl	64(%rsp), %eax
	movslq	24(%rsp), %r14
	movl	%eax, 20(%rsp)
	xorl	%eax, %eax
	leaq	(%r14,%r14,2), %r15
	salq	$2, %r14
	salq	$3, %r15
	.p2align 4,,10
	.p2align 3
.L12:
	movl	(%rcx,%r14), %edi
	movq	%r15, %rdx
	movq	%r15, %rsi
	addq	16(%rbx), %rdx
	movq	%rax, 8(%rsp)
	addq	8(%rbx), %rsi
	subl	$1, %edi
	call	_ZL11gpuirr_firriPdS_
	movq	(%rbx), %rcx
	movq	8(%rsp), %rax
	addq	$24, %r15
	movslq	(%rcx,%r14), %rdx
	addq	$4, %r14
	addq	$1, %rax
	subq	$1, %rdx
	imulq	$10016, %rdx, %rdx
	addl	12(%r13,%rdx), %r12d
	movl	24(%rsp), %edx
	addl	%eax, %edx
	cmpl	%edx, 20(%rsp)
	jg	.L12
	leaq	64(%rsp), %rsi
	leaq	32(%rsp), %rdi
	call	GOMP_loop_guided_next
	testb	%al, %al
	jne	.L13
.L11:
	call	GOMP_loop_end_nowait
	lock addl	%r12d, 28(%rbx)
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE1268:
	.size	_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.0, .-_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.0
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC3:
	.string	"gpuirr: it is already open\n"
.LC4:
	.string	"irrlib/gpuirr.avx.cpp"
.LC5:
	.string	"lmax <= 1 + NBlist::NB_MAX"
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC6:
	.string	"**************************** \n"
	.align 8
.LC7:
	.string	"Opening GPUIRR lib. AVX ver. \n"
	.section	.rodata.str1.1
.LC8:
	.string	" nmax = %d, lmax = %d\n"
	.section	.rodata.str1.8
	.align 8
.LC9:
	.string	"0 == posix_memalign(&ptr, 64, (1+nmax) * sizeof(Particle))"
	.align 8
.LC11:
	.string	"0 == posix_memalign(&ptr, 64, nmax * sizeof(NBlist))"
	.text
	.p2align 4,,15
	.globl	gpuirr_open_
	.type	gpuirr_open_, @function
gpuirr_open_:
.LFB1259:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$24, %rsp
	.cfi_def_cfa_offset 48
	cmpb	$0, _ZL7is_open(%rip)
	movl	(%rsi), %ebp
	movl	(%rdi), %ebx
	jne	.L23
	cmpl	$2501, %ebp
	jg	.L24
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC6, %edi
	movq	stderr(%rip), %rcx
	call	fwrite
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC7, %edi
	call	fwrite
	movl	%ebp, %r8d
	movq	stderr(%rip), %rdi
	movl	%ebx, %ecx
	movl	$.LC8, %edx
	movl	$1, %esi
	xorl	%eax, %eax
	call	__fprintf_chk
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC6, %edi
	call	fwrite
	leal	1(%rbx), %eax
	leaq	8(%rsp), %rdi
	cltq
	movl	$64, %esi
	leaq	(%rax,%rax,2), %rbp
	salq	$5, %rbp
	movq	%rbp, %rdx
	call	posix_memalign
	testl	%eax, %eax
	jne	.L25
	movq	8(%rsp), %rdi
	movq	%rbp, %rdx
	movl	$255, %esi
	call	memset
	movslq	%ebx, %rbp
	leaq	0(%rbp,%rbp,2), %rax
	movq	8(%rsp), %rdx
	movl	$64, %esi
	imulq	$10016, %rbp, %rbp
	vmovaps	.LC10(%rip), %xmm0
	salq	$5, %rax
	leaq	8(%rsp), %rdi
	addq	%rdx, %rax
	movq	%rdx, _ZL4ptcl(%rip)
	vmovaps	%xmm0, (%rax)
	vxorps	%xmm0, %xmm0, %xmm0
	movq	%rbp, %rdx
	vmovaps	%xmm0, 16(%rax)
	vmovaps	%xmm0, 32(%rax)
	vmovaps	%xmm0, 48(%rax)
	vmovaps	%xmm0, 64(%rax)
	vxorpd	%xmm0, %xmm0, %xmm0
	vmovapd	%xmm0, 80(%rax)
	call	posix_memalign
	testl	%eax, %eax
	jne	.L26
	movq	8(%rsp), %rdi
	movq	%rbp, %rdx
	movl	$255, %esi
	call	memset
	movq	8(%rsp), %rax
	xorl	%edx, %edx
	xorl	%esi, %esi
	movl	$_ZL11gpuirr_openii._omp_fn.1, %edi
	movl	%ebx, _ZL4nmax(%rip)
	movq	%rax, _ZL4list(%rip)
	call	GOMP_parallel_start
	xorl	%edi, %edi
	call	_ZL11gpuirr_openii._omp_fn.1
	call	GOMP_parallel_end
	movb	$1, _ZL7is_open(%rip)
	movq	$0, _ZL9time_grav(%rip)
	movq	$0, _ZL9num_steps(%rip)
	movq	$0, _ZL9num_fcall(%rip)
	movq	$0, _ZL9num_inter(%rip)
	addq	$24, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L23:
	.cfi_restore_state
	movl	$27, %edx
	movl	$1, %esi
	movl	$.LC3, %edi
	movq	stderr(%rip), %rcx
	call	fwrite
	addq	$24, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
.L26:
	.cfi_restore_state
	movl	$_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, %ecx
	movl	$376, %edx
	movl	$.LC4, %esi
	movl	$.LC11, %edi
	call	__assert_fail
.L25:
	movl	$_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, %ecx
	movl	$368, %edx
	movl	$.LC4, %esi
	movl	$.LC9, %edi
	call	__assert_fail
.L24:
	movl	$_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, %ecx
	movl	$358, %edx
	movl	$.LC4, %esi
	movl	$.LC5, %edi
	call	__assert_fail
	.cfi_endproc
.LFE1259:
	.size	gpuirr_open_, .-gpuirr_open_
	.section	.rodata.str1.1
.LC13:
	.string	"gpuirr: it is already close\n"
	.section	.rodata.str1.8
	.align 8
.LC15:
	.string	"Closing GPUIRR lib. CPU ver. \n"
	.section	.rodata.str1.1
.LC16:
	.string	"time grav  : %f sec\n"
.LC19:
	.string	"perf grav  : %f Gflops\n"
.LC20:
	.string	"perf grav  : %f usec\n"
.LC21:
	.string	"<#NB>      : %f \n"
	.text
	.p2align 4,,15
	.globl	gpuirr_close_
	.type	gpuirr_close_, @function
gpuirr_close_:
.LFB1260:
	.cfi_startproc
	subq	$136, %rsp
	.cfi_def_cfa_offset 144
	cmpb	$0, _ZL7is_open(%rip)
	je	.L35
	movq	_ZL4ptcl(%rip), %rdi
	call	free
	movq	_ZL4list(%rip), %rdi
	movq	$0, _ZL4ptcl(%rip)
	call	free
	movq	_ZL9num_inter(%rip), %rax
	movq	$0, _ZL4list(%rip)
	testq	%rax, %rax
	js	.L29
	vcvtsi2sdq	%rax, %xmm2, %xmm2
.L30:
	movq	_ZL9num_fcall(%rip), %rax
	vmovsd	_ZL9time_grav(%rip), %xmm1
	testq	%rax, %rax
	js	.L31
	vcvtsi2sdq	%rax, %xmm4, %xmm4
.L32:
	vdivsd	%xmm4, %xmm1, %xmm4
	movq	_ZL9num_steps(%rip), %rax
	testq	%rax, %rax
	vmulsd	.LC14(%rip), %xmm4, %xmm4
	js	.L33
	vcvtsi2sdq	%rax, %xmm3, %xmm3
.L34:
	vdivsd	%xmm3, %xmm2, %xmm3
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC6, %edi
	vmovsd	%xmm1, 64(%rsp)
	vmovsd	%xmm2, 96(%rsp)
	vmovsd	%xmm4, 32(%rsp)
	vmovsd	%xmm3, (%rsp)
	call	fwrite
	movl	$30, %edx
	movl	$1, %esi
	movq	stderr(%rip), %rcx
	movl	$.LC15, %edi
	call	fwrite
	movq	stderr(%rip), %rdi
	movl	$.LC16, %edx
	movl	$1, %esi
	movl	$1, %eax
	vmovsd	_ZL9time_grav(%rip), %xmm0
	call	__fprintf_chk
	movq	stderr(%rip), %rsi
	movl	$10, %edi
	call	fputc
	vmovsd	96(%rsp), %xmm2
	vmovsd	64(%rsp), %xmm1
	movl	$.LC19, %edx
	movl	$1, %esi
	vmulsd	.LC17(%rip), %xmm2, %xmm0
	movl	$1, %eax
	movq	stderr(%rip), %rdi
	vmulsd	.LC18(%rip), %xmm0, %xmm0
	vdivsd	%xmm1, %xmm0, %xmm0
	call	__fprintf_chk
	movq	stderr(%rip), %rdi
	vmovsd	32(%rsp), %xmm4
	movl	$.LC20, %edx
	movl	$1, %esi
	movl	$1, %eax
	vmovapd	%xmm4, %xmm0
	call	__fprintf_chk
	movq	stderr(%rip), %rdi
	vmovsd	(%rsp), %xmm3
	movl	$.LC21, %edx
	movl	$1, %esi
	movl	$1, %eax
	vmovapd	%xmm3, %xmm0
	call	__fprintf_chk
	movq	stderr(%rip), %rcx
	movl	$30, %edx
	movl	$1, %esi
	movl	$.LC6, %edi
	call	fwrite
	movb	$0, _ZL7is_open(%rip)
	addq	$136, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L35:
	.cfi_restore_state
	movl	$28, %edx
	movl	$1, %esi
	movl	$.LC13, %edi
	movq	stderr(%rip), %rcx
	addq	$136, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	jmp	fwrite
	.p2align 4,,10
	.p2align 3
.L29:
	.cfi_restore_state
	movq	%rax, %rdx
	andl	$1, %eax
	shrq	%rdx
	orq	%rax, %rdx
	vcvtsi2sdq	%rdx, %xmm2, %xmm2
	vaddsd	%xmm2, %xmm2, %xmm2
	jmp	.L30
	.p2align 4,,10
	.p2align 3
.L33:
	movq	%rax, %rdx
	andl	$1, %eax
	shrq	%rdx
	orq	%rax, %rdx
	vcvtsi2sdq	%rdx, %xmm3, %xmm3
	vaddsd	%xmm3, %xmm3, %xmm3
	jmp	.L34
	.p2align 4,,10
	.p2align 3
.L31:
	movq	%rax, %rdx
	andl	$1, %eax
	shrq	%rdx
	orq	%rax, %rdx
	vcvtsi2sdq	%rdx, %xmm4, %xmm4
	vaddsd	%xmm4, %xmm4, %xmm4
	jmp	.L32
	.cfi_endproc
.LFE1260:
	.size	gpuirr_close_, .-gpuirr_close_
	.p2align 4,,15
	.globl	gpuirr_set_jp_
	.type	gpuirr_set_jp_, @function
gpuirr_set_jp_:
.LFB1261:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-32, %rsp
	addq	$16, %rsp
	movl	(%rdi), %eax
	vmovsd	16(%rsi), %xmm1
	vmovsd	(%rsi), %xmm2
	vmovhpd	(%r9), %xmm1, %xmm0
	vmovsd	(%rdx), %xmm4
	vmovhpd	8(%rsi), %xmm2, %xmm1
	subl	$1, %eax
	vmovsd	(%rcx), %xmm5
	vinsertf128	$0x1, %xmm0, %ymm1, %ymm1
	cltq
	vcvtpd2psy	%ymm1, %xmm0
	vcvtps2pd	%xmm0, %ymm2
	vsubpd	%ymm2, %ymm1, %ymm1
	leaq	(%rax,%rax,2), %rax
	vmovsd	16(%rdx), %xmm2
	vmovhpd	8(%rdx), %xmm4, %xmm3
	salq	$5, %rax
	vmovhpd	8(%rcx), %xmm5, %xmm4
	vmovsd	(%r8), %xmm6
	vmovaps	%xmm0, %xmm0
	vinsertf128	$0x1, %xmm2, %ymm3, %ymm3
	movq	16(%rbp), %rdx
	vcvtpd2psy	%ymm3, %xmm3
	vmovsd	16(%rcx), %xmm2
	vmovhpd	8(%r8), %xmm6, %xmm5
	vmovaps	%xmm3, %xmm3
	addq	_ZL4ptcl(%rip), %rax
	vcvtpd2psy	%ymm1, %xmm1
	vinsertf128	$0x1, %xmm1, %ymm0, %ymm0
	vinsertf128	$0x1, %xmm2, %ymm4, %ymm4
	vmovsd	16(%r8), %xmm2
	vcvtpd2psy	%ymm4, %xmm4
	vmovddup	(%rdx), %xmm1
	vinsertf128	$0x1, %xmm4, %ymm3, %ymm3
	vinsertf128	$0x1, %xmm2, %ymm5, %ymm2
	vcvtpd2psy	%ymm2, %xmm2
	vmovaps	%xmm2, %xmm2
	vmovaps	%ymm0, (%rax)
	vinsertf128	$0x1, %xmm1, %ymm2, %ymm2
	vmovaps	%ymm3, 32(%rax)
	vmovaps	%ymm2, 64(%rax)
	leave
	.cfi_def_cfa 7, 8
	vzeroupper
	ret
	.cfi_endproc
.LFE1261:
	.size	gpuirr_set_jp_, .-gpuirr_set_jp_
	.section	.rodata.str1.1
.LC22:
	.string	"nnb <= NBlist::NB_MAX"
	.text
	.p2align 4,,15
	.globl	gpuirr_set_list_
	.type	gpuirr_set_list_, @function
gpuirr_set_list_:
.LFB1262:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	leaq	4(%rsi), %rcx
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movl	(%rsi), %eax
	movl	(%rdi), %r11d
	cmpl	$2500, %eax
	jg	.L51
	leal	-1(%r11), %r9d
	movq	_ZL4list(%rip), %r8
	movslq	%r9d, %r9
	imulq	$10016, %r9, %r9
	addq	%r8, %r9
	testl	%eax, %eax
	movl	%eax, 12(%r9)
	jle	.L39
	vmovdqa	.LC23(%rip), %xmm0
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L40:
	vmovdqu	(%rcx), %xmm2
	vmovdqu	16(%rcx), %xmm1
	addq	$32, %rcx
	vpsubd	%xmm0, %xmm2, %xmm2
	vpsubd	%xmm0, %xmm1, %xmm1
	vmovntdq	%xmm2, 16(%r9,%rdx,4)
	vmovntdq	%xmm1, 32(%r9,%rdx,4)
	addq	$8, %rdx
	cmpl	%edx, %eax
	jg	.L40
.L39:
	leal	6(%rax), %edi
	movl	%eax, %edx
	movl	_ZL4nmax(%rip), %esi
	subl	$1, %edx
	cmovns	%edx, %edi
	sarl	$3, %edi
	leal	8(,%rdi,8), %edi
	cmpl	%edi, %eax
	jge	.L37
	movslq	%r11d, %r11
	movl	%edi, %r10d
	imulq	$10016, %r11, %rbx
	subl	%eax, %r10d
	movl	%r10d, %ecx
	addq	%r8, %rbx
	movslq	%eax, %r8
	leaq	-2500(%r8), %rbp
	leaq	(%rbx,%rbp,4), %rdx
	andl	$15, %edx
	shrq	$2, %rdx
	negq	%rdx
	andl	$3, %edx
	cmpl	%r10d, %edx
	cmovbe	%edx, %ecx
	testl	%ecx, %ecx
	movl	%ecx, %r12d
	je	.L42
	leaq	(%r9,%r8,4), %r8
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L43:
	movl	%esi, 16(%r8,%rdx,4)
	addq	$1, %rdx
	addl	$1, %eax
	cmpl	%edx, %ecx
	ja	.L43
	cmpl	%ecx, %r10d
	je	.L37
.L42:
	subl	%ecx, %r10d
	movl	%r10d, %r8d
	shrl	$2, %r8d
	leal	0(,%r8,4), %r11d
	testl	%r11d, %r11d
	je	.L44
	addq	%r12, %rbp
	vmovd	%esi, %xmm1
	xorl	%edx, %edx
	leaq	(%rbx,%rbp,4), %rcx
	vpshufd	$0, %xmm1, %xmm0
	.p2align 4,,10
	.p2align 3
.L45:
	addl	$1, %edx
	vmovdqa	%xmm0, (%rcx)
	addq	$16, %rcx
	cmpl	%r8d, %edx
	jb	.L45
	addl	%r11d, %eax
	cmpl	%r11d, %r10d
	je	.L37
.L44:
	movslq	%eax, %rdx
	leaq	16(%r9,%rdx,4), %rdx
	.p2align 4,,10
	.p2align 3
.L46:
	addl	$1, %eax
	movl	%esi, (%rdx)
	addq	$4, %rdx
	cmpl	%eax, %edi
	jg	.L46
.L37:
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
.L51:
	.cfi_restore_state
	movl	$_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__, %ecx
	movl	$433, %edx
	movl	$.LC4, %esi
	movl	$.LC22, %edi
	call	__assert_fail
	.cfi_endproc
.LFE1262:
	.size	gpuirr_set_list_, .-gpuirr_set_list_
	.p2align 4,,15
	.globl	gpuirr_pred_all_
	.type	gpuirr_pred_all_, @function
gpuirr_pred_all_:
.LFB1263:
	.cfi_startproc
	vmovddup	(%rdx), %xmm0
	vmovapd	%xmm0, _ZL8vec_tnow(%rip)
	ret
	.cfi_endproc
.LFE1263:
	.size	gpuirr_pred_all_, .-gpuirr_pred_all_
	.p2align 4,,15
	.globl	gpuirr_pred_act_
	.type	gpuirr_pred_act_, @function
gpuirr_pred_act_:
.LFB1264:
	.cfi_startproc
	vmovddup	(%rdx), %xmm0
	vmovapd	%xmm0, _ZL8vec_tnow(%rip)
	ret
	.cfi_endproc
.LFE1264:
	.size	gpuirr_pred_act_, .-gpuirr_pred_act_
	.p2align 4,,15
	.globl	gpuirr_firr_vec_
	.type	gpuirr_firr_vec_, @function
gpuirr_firr_vec_:
.LFB1265:
	.cfi_startproc
	movq	%rbx, -32(%rsp)
	movq	%rsi, %rbx
	.cfi_offset 3, -40
	xorl	%esi, %esi
	movq	%rbp, -24(%rsp)
	movq	%r12, -16(%rsp)
	movq	%rdx, %r12
	.cfi_offset 12, -24
	.cfi_offset 6, -32
	movq	%r13, -8(%rsp)
	subq	$136, %rsp
	.cfi_def_cfa_offset 144
	.cfi_offset 13, -16
	movq	%rcx, %r13
	movslq	(%rdi), %rbp
	leaq	80(%rsp), %rdi
	call	gettimeofday
	xorl	%edx, %edx
	vcvtsi2sdq	88(%rsp), %xmm1, %xmm1
	leaq	48(%rsp), %rsi
	movl	$_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.0, %edi
	vmovsd	.LC24(%rip), %xmm0
	movl	$0, 76(%rsp)
	vcvtsi2sdq	80(%rsp), %xmm2, %xmm2
	movl	%ebp, 72(%rsp)
	movq	%rbx, 48(%rsp)
	vmulsd	%xmm0, %xmm1, %xmm1
	movq	%r12, 56(%rsp)
	movq	%r13, 64(%rsp)
	vmovsd	%xmm0, (%rsp)
	vaddsd	%xmm1, %xmm2, %xmm1
	vmovsd	%xmm1, 40(%rsp)
	call	GOMP_parallel_start
	leaq	48(%rsp), %rdi
	call	_ZL15gpuirr_firr_veciPKiPA3_dS2_._omp_fn.0
	call	GOMP_parallel_end
	movslq	76(%rsp), %rax
	leaq	80(%rsp), %rdi
	xorl	%esi, %esi
	addq	%rax, _ZL9num_inter(%rip)
	call	gettimeofday
	vmovsd	(%rsp), %xmm0
	vcvtsi2sdq	88(%rsp), %xmm2, %xmm2
	vcvtsi2sdq	80(%rsp), %xmm1, %xmm1
	addq	%rbp, _ZL9num_steps(%rip)
	addq	$1, _ZL9num_fcall(%rip)
	movq	104(%rsp), %rbx
	vmulsd	%xmm0, %xmm2, %xmm0
	movq	112(%rsp), %rbp
	movq	120(%rsp), %r12
	movq	128(%rsp), %r13
	vaddsd	%xmm0, %xmm1, %xmm0
	vsubsd	40(%rsp), %xmm0, %xmm0
	vaddsd	_ZL9time_grav(%rip), %xmm0, %xmm0
	vmovsd	%xmm0, _ZL9time_grav(%rip)
	addq	$136, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1265:
	.size	gpuirr_firr_vec_, .-gpuirr_firr_vec_
	.section	.text.startup,"ax",@progbits
	.p2align 4,,15
	.type	_GLOBAL__sub_I_gpuirr_open_, @function
_GLOBAL__sub_I_gpuirr_open_:
.LFB1267:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	movl	$_ZStL8__ioinit, %edi
	call	_ZNSt8ios_base4InitC1Ev
	movl	$__dso_handle, %edx
	movl	$_ZStL8__ioinit, %esi
	movl	$_ZNSt8ios_base4InitD1Ev, %edi
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	jmp	__cxa_atexit
	.cfi_endproc
.LFE1267:
	.size	_GLOBAL__sub_I_gpuirr_open_, .-_GLOBAL__sub_I_gpuirr_open_
	.section	.ctors,"aw",@progbits
	.align 8
	.quad	_GLOBAL__sub_I_gpuirr_open_
	.local	_ZL7is_open
	.comm	_ZL7is_open,1,1
	.section	.rodata
	.align 16
	.type	_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, @object
	.size	_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__, 27
_ZZL11gpuirr_openiiE19__PRETTY_FUNCTION__:
	.string	"void gpuirr_open(int, int)"
	.local	_ZL4ptcl
	.comm	_ZL4ptcl,8,32
	.local	_ZL4list
	.comm	_ZL4list,8,32
	.local	_ZL4nmax
	.comm	_ZL4nmax,4,32
	.local	_ZL9time_grav
	.comm	_ZL9time_grav,8,8
	.local	_ZL9num_steps
	.comm	_ZL9num_steps,8,8
	.local	_ZL9num_fcall
	.comm	_ZL9num_fcall,8,8
	.local	_ZL9num_inter
	.comm	_ZL9num_inter,8,8
	.align 32
	.type	_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__, @object
	.size	_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__, 43
_ZZL15gpuirr_set_listiiPKiE19__PRETTY_FUNCTION__:
	.string	"void gpuirr_set_list(int, int, const int*)"
	.local	_ZL8vec_tnow
	.comm	_ZL8vec_tnow,16,16
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.local	_ZL11num_threads
	.comm	_ZL11num_threads,4,32
	.weakref	_ZL20__gthrw_pthread_oncePiPFvvE,pthread_once
	.weakref	_ZL27__gthrw_pthread_getspecificj,pthread_getspecific
	.weakref	_ZL27__gthrw_pthread_setspecificjPKv,pthread_setspecific
	.weakref	_ZL22__gthrw_pthread_createPmPK14pthread_attr_tPFPvS3_ES3_,pthread_create
	.weakref	_ZL20__gthrw_pthread_joinmPPv,pthread_join
	.weakref	_ZL21__gthrw_pthread_equalmm,pthread_equal
	.weakref	_ZL20__gthrw_pthread_selfv,pthread_self
	.weakref	_ZL22__gthrw_pthread_detachm,pthread_detach
	.weakref	_ZL22__gthrw_pthread_cancelm,pthread_cancel
	.weakref	_ZL19__gthrw_sched_yieldv,sched_yield
	.weakref	_ZL26__gthrw_pthread_mutex_lockP15pthread_mutex_t,pthread_mutex_lock
	.weakref	_ZL29__gthrw_pthread_mutex_trylockP15pthread_mutex_t,pthread_mutex_trylock
	.weakref	_ZL31__gthrw_pthread_mutex_timedlockP15pthread_mutex_tPK8timespec,pthread_mutex_timedlock
	.weakref	_ZL28__gthrw_pthread_mutex_unlockP15pthread_mutex_t,pthread_mutex_unlock
	.weakref	_ZL26__gthrw_pthread_mutex_initP15pthread_mutex_tPK19pthread_mutexattr_t,pthread_mutex_init
	.weakref	_ZL29__gthrw_pthread_mutex_destroyP15pthread_mutex_t,pthread_mutex_destroy
	.weakref	_ZL30__gthrw_pthread_cond_broadcastP14pthread_cond_t,pthread_cond_broadcast
	.weakref	_ZL27__gthrw_pthread_cond_signalP14pthread_cond_t,pthread_cond_signal
	.weakref	_ZL25__gthrw_pthread_cond_waitP14pthread_cond_tP15pthread_mutex_t,pthread_cond_wait
	.weakref	_ZL30__gthrw_pthread_cond_timedwaitP14pthread_cond_tP15pthread_mutex_tPK8timespec,pthread_cond_timedwait
	.weakref	_ZL28__gthrw_pthread_cond_destroyP14pthread_cond_t,pthread_cond_destroy
	.weakref	_ZL26__gthrw_pthread_key_createPjPFvPvE,pthread_key_create
	.weakref	_ZL26__gthrw_pthread_key_deletej,pthread_key_delete
	.weakref	_ZL30__gthrw_pthread_mutexattr_initP19pthread_mutexattr_t,pthread_mutexattr_init
	.weakref	_ZL33__gthrw_pthread_mutexattr_settypeP19pthread_mutexattr_ti,pthread_mutexattr_settype
	.weakref	_ZL33__gthrw_pthread_mutexattr_destroyP19pthread_mutexattr_t,pthread_mutexattr_destroy
	.section	.rodata.cst32,"aM",@progbits,32
	.align 32
.LC0:
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.long	1069547520
	.align 32
.LC1:
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.long	3204448256
	.align 32
.LC2:
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.long	3225419776
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC10:
	.long	1132396544
	.long	1132396544
	.long	1132396544
	.long	0
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC14:
	.long	0
	.long	1093567616
	.align 8
.LC17:
	.long	0
	.long	1078853632
	.align 8
.LC18:
	.long	3894859413
	.long	1041313291
	.section	.rodata.cst16
	.align 16
.LC23:
	.long	1
	.long	1
	.long	1
	.long	1
	.section	.rodata.cst8
	.align 8
.LC24:
	.long	2696277389
	.long	1051772663
	.ident	"GCC: (Ubuntu/Linaro 4.6.3-1ubuntu5) 4.6.3"
	.section	.note.GNU-stack,"",@progbits
