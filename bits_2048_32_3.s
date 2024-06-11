	.file	"bits_2048_32_3.c"
    .section        .rodata
    .align 16
.LC0:
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
    .byte   1
.LC1:
    .byte   255
    .byte   0
    .byte   255
    .byte   0
    .byte   255
    .byte   0
    .byte   255
    .byte   0
    .byte   255
    .byte   0
    .byte   255
    .byte   0
    .byte   255
    .byte   0
    .byte   255
    .byte   0
	.text
.globl bits_2048_32_3
	.type	bits_2048_32_3, @function
bits_2048_32_3:
	pushl	%ebp
    movl	%esp, %ebp
	movl	8(%ebp), %eax #one
	movl	12(%ebp), %eax #two
	movl	16(%ebp), %eax #thr
	movl	32(%ebp), %eax #rem
    pxor    %xmm7, %xmm7       #   Clear xmm7
    pxor    %xmm6, %xmm6
	movdqa  .LC0, %xmm0        #   Load the 2-bit bitmask
	movl    32(%ebp), %eax     #   Pointer to remainder
	movq    (%eax), %xmm1
    movl    36(%ebp), %eax     #   Pointer to data[0]
	lddqu   (%eax), %xmm2      #   Read in data
	punpcklqdq   %xmm2, %xmm1  #   Unite the remainder of prev. data with start of new data
    movl	$0, %eax
	leal	(%eax,%eax), %edx  #   set the offset for destination arrays
	movdqa  %xmm1, %xmm2   #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2  #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2   #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2   #   clean up roll-overs from subtraction
	movdqa  %xmm2, %xmm3   #   store a copy of one[i] data in xmm3
	psrlw   $2, %xmm1      #   Left shift the data by two bits
    movdqa  %xmm1, %xmm2   #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2      #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2   #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2   #   clean up roll-overs from subtraction
	movdqa  %xmm2, %xmm4   #   store a copy of two[i] data in xmm4
	psrlw   $2, %xmm1      #   Left shift the data by two bits
    movdqa  %xmm1, %xmm2   #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2      #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2   #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2   #   clean up roll-overs from subtraction
	movdqa  %xmm2, %xmm5   #   store a copy of one[i] data in xmm5

    movl	$1, %eax
	leal	(%eax,%eax), %edx #   set the offset for destination arrays
    movl	36(%ebp), %eax  #   Pointer to data[0]
	lddqu   (%eax), %xmm1     #    Load counter data to xmm1
	movdqa  %xmm1, %xmm2    #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2       #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2    #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2    #   clean up roll-overs from subtraction
    punpcklbw %xmm3, %xmm2  #   combine 8+8 pieces of output for summing
    psadbw  %xmm7, %xmm2    # sum 16xbytes into two words
    movdqa  %xmm2, %xmm3    # sum the two words
    psrldq  $8, %xmm3       # sum the two words
    paddq   %xmm3, %xmm2    # sum the two words
    pextrw  $0, %xmm2, %edx # extract the sum into a 32bit-register
    movl	8(%ebp), %eax   #   Get one[i] array address
    movw    %dx, (%eax)     # copy the sum into memory (short)
	psrlw   $2, %xmm1       #   Left shift the data by two bits
    movdqa  %xmm1, %xmm2    #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2       #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2    #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2    #   clean up roll-overs from subtraction
    punpcklbw %xmm4, %xmm2  #   combine 8+8 pieces of output for summing
    psadbw  %xmm7, %xmm2    # sum 16xbytes into two words
    movdqa  %xmm2, %xmm4    # sum the two words
    psrldq  $8, %xmm4       # sum the two words
    paddq   %xmm4, %xmm2    # sum the two words
    pextrw  $0, %xmm2, %edx # extract the sum into a 32bit-register
    movl	12(%ebp), %eax  #   Get two[i] array address
    movw    %dx, (%eax)     # copy the sum into memory (short)
	psrlw   $2, %xmm1       #   Left shift the data by two bits
    movdqa  %xmm1, %xmm2    #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2       #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2    #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2    #   clean up roll-overs from subtraction
	punpcklbw %xmm5, %xmm2  #   combine 8+8 pieces of output for summing
    psadbw  %xmm7, %xmm2    # sum 16xbytes into two words
    movdqa  %xmm2, %xmm5    # sum the two words
    psrldq  $8, %xmm5       # sum the two words
    paddq   %xmm5, %xmm2    # sum the two words
    pextrw  $0, %xmm2, %edx  # extract the sum into a 32bit-register
    movl	16(%ebp), %eax     #   Get thr[i] array address
    movw    %dx, (%eax)     # copy the sum into memory (short)
    movl    $1, -8(%ebp)    # one[i] i value
	movl	$1, -4(%ebp)    # data[i] i value

.L2:
	movl	-4(%ebp), %eax
	cmpl	$255, %eax
	jl	.L5
	jmp	.L3
.L5:
	leal	(%eax,%eax), %edx #   set the offset for data array
	movl	36(%ebp), %eax    #   Get counter data address
	lddqu   (%eax,%edx,4), %xmm1     #    Load counter data to xmm1
    leal	-4(%ebp), %eax     #   Increment counter
	incl	(%eax)
	movdqa  %xmm1, %xmm2   #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2      #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2   #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2   #   clean up roll-overs from subtraction
	movdqa  %xmm2, %xmm3   #   store a copy of one[i] data in xmm3
	psrlw   $2, %xmm1      #   Left shift the data by two bits
    movdqa  %xmm1, %xmm2   #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2      #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2   #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2   #   clean up roll-overs from subtraction
	movdqa  %xmm2, %xmm4   #   store a copy of two[i] data in xmm4
	psrlw   $2, %xmm1      #   Left shift the data by two bits
    movdqa  %xmm1, %xmm2   #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2      #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2   #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2   #   clean up roll-overs from subtraction
	movdqa  %xmm2, %xmm5   #   store a copy of one[i] data in xmm5
	movl	-4(%ebp), %eax
	leal	(%eax,%eax), %edx #   set the offset for data array
	movl	36(%ebp), %eax    #   Get counter data address
	lddqu   (%eax,%edx,4), %xmm1     #    Load counter data to xmm1
    leal	-4(%ebp), %eax     #   Increment counter
	incl	(%eax)
	movdqa  %xmm1, %xmm2   #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2  #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2   #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2   #   clean up roll-overs from subtraction
    punpcklbw %xmm3, %xmm2 #   combine 8+8 pieces of output for summing
    psadbw  %xmm7, %xmm2    # sum 16xbytes into two words
    movdqa  %xmm2, %xmm3    # sum the two words
    psrldq  $8, %xmm3       # sum the two words
    paddq   %xmm3, %xmm2    # sum the two words
    pextrw  $0, %xmm2, %ecx  # extract the sum into a 32bit-register
	movl	-8(%ebp), %edx # leal	(%eax,%eax), %edx #   set the offset for data array
    leal	-8(%ebp), %eax     #   Increment counter
	incl	(%eax)
    movl	8(%ebp), %eax     #   Get one[i] array address
    movb    %cl, (%eax,%edx,1)     # copy the sum into memory (short)
	psrlw   $2, %xmm1      #   Left shift the data by two bits
    movdqa  %xmm1, %xmm2   #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2  #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2   #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2   #   clean up roll-overs from subtraction
    punpcklbw %xmm4, %xmm2 #   combine 8+8 pieces of output for summing
    psadbw  %xmm7, %xmm2    # sum 16xbytes into two words
    movdqa  %xmm2, %xmm4    # sum the two words
    psrldq  $8, %xmm4       # sum the two words
    paddq   %xmm4, %xmm2    # sum the two words
    pextrw  $0, %xmm2, %ecx  # extract the sum into a 32bit-register
    movl	12(%ebp), %eax     #   Get two[i] array address
    movb    %cl, (%eax,%edx,1)     # copy the sum into memory (short)
	psrlw   $2, %xmm1      #   Left shift the data by two bits
    movdqa  %xmm1, %xmm2   #   Create a copy of counter data in xmm2
	psrldq  $1, %xmm2  #   shift xmm2 by one byte (time + 1)
	psubb   %xmm1, %xmm2   #   counter(t+1) - counter(t)
	pand    %xmm0, %xmm2   #   clean up roll-overs from subtraction
	punpcklbw %xmm5, %xmm2 #   combine 8+8 pieces of output for summing
    psadbw  %xmm7, %xmm2    # sum 16xbytes into two words
    movdqa  %xmm2, %xmm5    # sum the two words
    psrldq  $8, %xmm5       # sum the two words
    paddq   %xmm5, %xmm2    # sum the two words
    pextrw  $0, %xmm2, %ecx  # extract the sum into a 32bit-register
    movl	16(%ebp), %eax     #   Get thr[i] array address
    movb    %cl, (%eax,%edx,1)     # copy the sum into memory (short)
	jmp	.L2
.L3:

    pxor    %xmm6, %xmm6
    movdqa  .LC1, %xmm5
    
    movl	8(%ebp), %eax    #    Get scratch memory address for colour one
    movl    $0, %edx
    movl	20(%ebp), %ecx
    
 	lddqu   (%eax, %edx,8), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx,8), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx 	
    lddqu   (%eax, %edx,8), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx,8), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
    lddqu   (%eax, %edx,8), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx,8), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx 	
    lddqu   (%eax, %edx,8), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx,8), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)

    movl	12(%ebp), %eax    #    Get scratch memory address for colour two
    movl    $0, %edx
    movl	24(%ebp), %ecx
    
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx 	
    lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
    lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx 	
    lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)

    movl	16(%ebp), %eax    #    Get scratch memory address for colour thr
    movl    $0, %edx
    movl	28(%ebp), %ecx
    
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx 	
    lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
    lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx 	
    lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)
    incl    %edx
    incl    %edx
 	lddqu   (%eax, %edx), %xmm0     #    scratch data to xmm1
    movdqa  %xmm0, %xmm7
    psrldq  $1, %xmm7
    paddb   %xmm7, %xmm0
    pand    %xmm5, %xmm0
    movdqu  %xmm0, (%ecx, %edx,8)

    movl	$254, %eax
    leal	(%eax,%eax), %edx #   access the data[last-1]
    movl	36(%ebp), %eax    #   Get data address
	lddqu   (%eax, %edx, 4), %xmm1     #    Load data to xmm1
    psrldq  $8, %xmm1   #   shift xmm1 by 8 bytes
    movl    32(%ebp), %eax  #Add data to rem(ainder)
    movq    %xmm1, (%eax)
	popl	%ebp
	ret
	.size	bits_2048_32_3, .-bits_2048_32_3
	.section	.note.GNU-stack,"",@progbits
    .ident  "Frogcode"
