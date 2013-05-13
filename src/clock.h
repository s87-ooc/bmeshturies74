#ifndef __CLOCK_H__
#define __CLOCK_H__

/*************************************************************

 Projet

 (C) 2013 Charles Podkanski (charles@podkanski.com),
          Stjepan Stamenkovic (stjepan@stjepan.net)

 ---

     Fichier: clock.h

 Description: macros pour mesurer le temps du code

**************************************************************/

clock_t _t;

#define RESETCLOCK _t = clock
#define CLOCK(T) T = clock() - _t
#define LOGTIME(S, T)cout << S << ": " << ((double)T) / CLOCKS_PER_SEC << "s (" << T << " clocks)" << endl;
#define LOGPARTTIME(S, T, TMAX)cout << S << ": " << ((double)T) / CLOCKS_PER_SEC << "s (" << T << " clocks) "<< (double)T / TMAX * 100. << "%" << endl;

#endif // __CLOCK_H__
