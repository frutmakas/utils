/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/02/11 21:48:42 $
 * $Revision: 1.1.2.4 $
 * $Id: pilot.h,v 1.1.2.4 2004/02/11 21:48:42 syed Exp $
 ********************************************************************/
#ifndef _TOOLS_PILOT_H_
#define _TOOLS_PILOT_H_
#include "tools/utilitis.h" 

ZVector insert_pilot(const ZVector &pilot, const ZVector &data, int symbols_between_pilot); 
ZVector delete_pilot(const ZVector data, int pilot_size, int symbols_between_pilot);
ZVector extract_pilot(const ZVector data, int pilot_size, int symbols_between_pilot) ;
int separate_pilot(const ZVector data, int pilot_size, int symbols_between_pilot, ZVector &pilot, ZVector &raw_data) ;

#endif

