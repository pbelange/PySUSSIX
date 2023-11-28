
# !=======================================================================
# !
# !  ORDERS THE HARMONICS FOUND BY SPECTRUM
# !
# !  INPUT PARAMETERS:
# !
# !  NARM : NUMBER OF HARMONICS TO BE OREDERED
# !  NR   : MAXIMUM ORDER OF HARMONIC TO BE LOOKED FOR IN THE LINEAR
# !         COMBINATIONS. (NR<10 ALTRIMENTI L'OUTPUT PUO INCASINARSI)
# !  EPS  : MAXIMUM ERROR ACCEPTED TO FIND THE LINEAR COMBINATIONS.
# !  TUNEX,TUNEZ,TUNES ARE THE EXPECTED TUNES WITHIN 0.01
# !
# !  MIND :  THE INPUT DATA COME FROM THE COMMON, SO THIS ROUTINE
# !          MUST BE USED ONLY AFTER THE CALL TO DATSPE
# !          THE ARRAYS TX,TY,TZ ARE USED IN ORDER NOT TO CHANGE
# !          THE ARRAYS TXA,TYA,TSA.
# !
# !  THE OUTPUT IS PLACED IN THE FILE lines WHICH IS NOT CLOSED AT THE
# !  END OF THE SUBROUTINE IN ORDER TO COLLECT THE RESULTS OF DIFFERENT
# !  SIGNALS IN ONE FILE ONLY. THE DIFFERENT SIGNALS ARE NUMBERED BY IOU.
# !
# !  A CHECK IS PERFORMED TO SEE IF THE HARMONICS ARE CLOSER THAN 1/NTURN
# !
# !  IOU  : INDEX WHICH IDENTIFIES THE CASE ANALIZED
# !
# !  AUTHOR: R.BARTOLINI 12/02/1996
# !  LAST MODIFIED: 21/09/1997
# !
# !=======================================================================
# ! ==============================================================================
# ! 2012 KEVIN: CHANGES MADE:
# ! ADDED OUTPUT VARIABLES: AMPLITUDE, PHASE, OX, AX, OY, AY, OS, AS
# ! COMMENTED MOST WRITE(30)
# ! REPLACED ALL write(30,100) -> CORRESPONDING VARIABLE ASSIGNMENTS I.E. OX, AX
# ! INSERTED SOME LINES FROM ROGELIO'S SUSSIX4DRIVEXX
# ! ==============================================================================

def ordres(eps, narm, nr, idam, iunit, nturn, tunex,
           tuney, tunez, istune, etune,
           amplitude, phase, ox, ax, oy, ay, os, as):
    
    # INITIALIZATION
    pi = 4.0 * np.arctan(1.0)
    checkn = 1.0 / float(nturn)
    imissx = 0
    imissy = 0
    imissz = 0
    iscax = 0
    iscay = 0
    iscaz = 0
    
    if nr > 10:
        print('ERROR IN ORDRES: NR LARGER THAN 10')
        return
    
    tx = np.zeros(narm)
    ty = np.zeros(narm)
    tz = np.zeros(narm)
    
    for j in range(1, narm + 1):
        tx[j - 1] = txa[j - 1]
        ty[j - 1] = tya[j - 1]
        tz[j - 1] = tsa[j - 1]

    # TUNES PARAMETERS AND EVENTUAL CHECK FOR THE EXPECTED TUNES
    if istune >= 1:
        # Check X TUNE
        dtunex = np.abs(np.abs(tx[0]) - np.abs(tunex))
        if dtunex > etune[0] or (tx[0] * tunex) < 0:
            print('X TUNE DIFFERENT FROM EXPECTED')
            print(-tx[0], -tunex)
            ntx = 1
            for nt in range(2, narm + 1):
                dtunex = np.abs(np.abs(tx[nt - 1]) - np.abs(tunex))
                if dtunex <= etune[0] and (tx[nt - 1] * tunex) > 0:
                    ntx = nt
                    break
            if ntx > 1:
                print('EXPECTED TUNE X FOUND AT LINE', ntx)
            elif ntx == 1:
                print('EXPECTED TUNE X NOT FOUND')
                if istune == 1:
                    print('LINE 1 ASSUMED AS TUNE!!!')
            tx[0], tx[ntx - 1] = tx[ntx - 1], txa[0]
        txt = tunex if istune == 2 else tx[0]
        pxt = np.abs(zxpes[0])
        if np.abs(pxt) > 0:
            pxtr = np.real(zxpes[0]) / pxt
            pxti = np.imag(zxpes[0]) / pxt
            if pxti == 0 and pxtr == 0:
                fxt = 0
            else:
                fxt = np.arctan2(pxti, pxtr)
        else:
            pxtr = 0
            pxti = 0
            fxt = 0
        fxt = fxt / pi * 180.0
        tyt = 9999.0
        tzt = 8888.0
        if idam >= 2:
            if istune >= 1:
                # Check Y TUNE
                dtuney = np.abs(np.abs(ty[0]) - np.abs(tuney))
                if dtuney > etune[1] or (ty[0] * tuney) < 0:
                    print('Y TUNE DIFFERENT FROM EXPECTED')
                    print(-ty[0], -tuney)
                    nty = 1
                    for nt in range(2, narm + 1):
                        dtuney = np.abs(np.abs(ty[nt - 1]) - np.abs(tuney))
                        if dtuney <= etune[1] and (ty[nt - 1] * tuney) > 0:
                            nty = nt
                            break
                    if nty > 1:
                        print('EXPECTED TUNE Y FOUND AT LINE', nty)
                    elif nty == 1:
                        print('EXPECTED TUNE Y NOT FOUND')
                        if istune == 1:
                            print('LINE 1 ASSUMED AS TUNE!!!')
                    ty[0], ty[nty - 1] = ty[nty - 1], tya[0]
                if istune == 2:
                    tyt = tuney
                else:
                    tyt = ty[0]
                pyt = np.abs(zypes[0])
                if np.abs(pyt) > 0:
                    pytr = np.real(zypes[0]) / pyt
                    pyti = np.imag(zypes[0]) / pyt
                    if pyti == 0 and pytr == 0:
                        fyt = 0
                    else:
                        fyt = np.arctan2(pyti, pytr)
                else:
                    pytr = 0
                    pyti = 0
                    fyt = 0
                fyt = fyt / pi * 180.0
                dty = np.abs(tyt - txt)
                if dty <= eps:
                    print('TUNEX AND TUNEY ARE TOO CLOSE')
                    tyt = 9999.0
            elif idam < 2 and istune == 2:
                tyt = tuney
        if idam == 3:
            if istune >= 1:
                # Check Z TUNE
                dtunez = np.abs(np.abs(tz[0]) - np.abs(tunez))
                if dtunez > etune[2] or (tz[0] * tunez) < 0:
                    print('Y TUNE DIFFERENT FROM EXPECTED')
                    print(-tz[0], -tunez)
                    ntz = 1
                    for nt in range(2, narm + 1):
                        dtunez = np.abs(np.abs(tz[nt - 1]) - np.abs(tunez))
                        if dtunez <= etune[2] and (tz[nt - 1] * tunez) > 0:
                            ntz = nt
                            break
                    if ntz > 1:
                        print('EXPECTED TUNE S FOUND AT LINE', ntz)
                    elif ntz == 1:
                        print('EXPECTED TUNE S NOT FOUND')
                        if istune == 1:
                            print('LINE 1 ASSUMED AS TUNE!!!')
                    tz[0], tz[ntz - 1] = tz[ntz - 1], tsa[0]
                if istune == 2:
                    tzt = tunez
                else:
                    tzt = tz[0]
                pzt = np.abs(ztpes[0])
                if np.abs(pzt) > 0:
                    pztr = np.real(ztpes[0]) / pzt
                    pzti = np.imag(ztpes[0]) / pzt
                    if pzti == 0 and pztr == 0:
                        fzt = 0
                    else:
                        fzt = np.arctan2(pzti, pztr)
                else:
                    pztr = 0
                    pzti = 0
                    fzt = 0
                fzt = fzt / pi * 180.0
                dtz = np.abs(tzt - txt)
                if dtz <= eps:
                    print('TUNEX AND TUNES ARE TOO CLOSE')
                    tzt = 8888.0
            elif idam < 2 and istune == 2:
                tzt = tunez

    # DAMPING
    for i in range(1, narm + 1):
        pxt = zxpes[i - 1]
        pyt = zypes[i - 1]
        pzt = ztpes[i - 1]
        psx = np.abs(pxt)
        psy = np.abs(pyt)
        psz = np.abs(pzt)
        if psx > 0 and os[i - 1] > 0:
            if idam == 1 or idam == 3:
                gxr = -pxtr / os[i - 1]
                gxi = -pxti / os[i - 1]
                gx = gxr + 1j * gxi
                ax[i - 1] = np.abs(gx)
                ox[i - 1] = np.angle(gx)
                if ax[i - 1] > 1:
                    iscax = 1
            else:
                ax[i - 1] = 0
                ox[i - 1] = 0
                os[i - 1] = 0
        else:
            ax[i - 1] = 0
            ox[i - 1] = 0
            os[i - 1] = 0
        if psy > 0 and os[i - 1] > 0:
            if idam == 1 or idam == 3:
                gyr = -pytr / os[i - 1]
                gyi = -pyti / os[i - 1]
                gy = gyr + 1j * gyi
                ay[i - 1] = np.abs(gy)
                oy[i - 1] = np.angle(gy)
                if ay[i - 1] > 1:
                    iscay = 1
            else:
                ay[i - 1] = 0
                oy[i - 1] = 0
                os[i - 1] = 0
        else:
            ay[i - 1] = 0
            oy[i - 1] = 0
            os[i - 1] = 0
        if psz > 0 and os[i - 1] > 0:
            if idam == 1 or idam == 3:
                gzr = -pztr / os[i - 1]
                gzi = -pzti / os[i - 1]
                gz = gzr + 1j * gzi
                az[i - 1] = np.abs(gz)
                oz[i - 1] = np.angle(gz)
                if az[i - 1] > 1:
                    iscaz = 1
            else:
                az[i - 1] = 0
                oz[i - 1] = 0
                os[i - 1] = 0
        else:
            az[i - 1] = 0
            oz[i - 1] = 0
            os[i - 1] = 0

    if idam == 1 or idam == 3:
        iscat = iscax + iscay + iscaz
        if iscat > 0:
            print('MODULI OF A PARAMETERS TOO LARGE')
            print('ARM MODULI AND PHASES SET TO ZERO')
            for j in range(1, narm + 1):
                ax[j - 1] = 0
                ay[j - 1] = 0
                az[j - 1] = 0
                ox[j - 1] = 0
                oy[j - 1] = 0
                oz[j - 1] = 0

    if idam >= 2:
        print('PHASES AND MODULI OF DAMPING FACTORS')
        print('ARM            MODULI          PHASES')
        for i in range(1, narm + 1):
            print(f'{i}   {os[i - 1]}   {ax[i - 1]}   {ay[i - 1]}   {az[i - 1]}   {ox[i - 1]}   {oy[i - 1]}   {oz[i - 1]}')

    return
