/* Copyright (C) 2010, Rodrigo Cánovas, all rights reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */



#include	<sys/types.h>
#include	<sys/times.h>
#include	<sys/param.h>


static struct tms start3;

void startclock3() {
  times(&start3);
}

long stopclock3() {
  struct tms t;
  times(&t);
  return(t.tms_utime - start3.tms_utime);
}

