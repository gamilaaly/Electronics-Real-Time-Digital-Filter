#include <compiler_defs.h>
#include <C8051F020_defs.h>


void ADC0_function(void){
	REF0CN = 0x03;
	ADC0CF = 0x00;
	AMX0CF = 0x00; 
	AMX0SL = 0x02;
	ADC0CN = 0x81; 
}

void DAC0_function(void){
	DAC0CN = 0x84;
	
}

unsigned int input[2] = {0,0};


void init_intClock(void)
{
	OSCXCN = 0x00;
	OSCICN = 0x14; // should be put as 0x14 to make it work	 but originally is it 0x04
	while (!(OSCICN & 0x10)); // Wait till IFRDY ( bit no 4 ) pin is set	, not working ,why if ^4 is 0
	/* Missing Clock Enable Bit(7) is 0 ( disabled) 
	  bits [ 1 : 0] are 10 for 2MHz
	  bit (2) is 1 to enable internal clock
	  bit (3) is 0 to disable external clock
	  */
}


void main(void){
EA = 0; 
WDTCN = 0xde;
WDTCN = 0xad;
init_intClock();
ADC0_function();
DAC0_function();
EA = 1;
XBR2 = 0x40;  
P1MDOUT = 0xff;
P1 = 0x00;

	while(1){
	AD0BUSY = 1;
	while(!AD0INT);
	AD0INT = 0;
	input[0] = input[1];
	input[1]= ADC0;
	DAC0 = (input[1] + input[0]);

	}
}