#include <compiler_defs.h>
#include <C8051F020_defs.h> // Include SFR declarations


void init_ADC0(void){
	REF0CN = 0x03; // Vref...use Internal Vref (2.4v)
	ADC0CF = 0x00; // gain = 1 ... Conversion Clock same as System Clock
	AMX0CF = 0x00; //
	AMX0SL = 0x02; // take input from Channel 2
	ADC0CN = 0x81; // enable ADC ... start conversion when AD0BUSY=1
}

void init_DAC0(void){
	DAC0CN = 0x84;
	
}

unsigned int x[2] = {0,0};


void init_intClock(void)
{
	OSCXCN = 0x00;
	// Bits 6-4 is 000 to disable it , else > don't care
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
init_ADC0();
init_DAC0();
EA = 1;
XBR2 = 0x40;  
P1MDOUT = 0xff;
P1 = 0x00;

	while(1){
	AD0BUSY = 1;
	while(!AD0INT); // becomes 1 when conversion ends
	AD0INT = 0;
	x[0] = x[1];
	x[1]= ADC0;
	DAC0 = (x[1] + x[0]);

	}
}
