program TopologyParser
    implicit none
    
    integer :: npx, npy, npz, npxyz, natoms, attract_label, nbas, nmol, num, iostat
    integer :: jobType = 0, i = 1, j = 1, indexList = 1, selection = 1, newIndex = 666
    double precision xmin, xmax, ymin, ymax, zmin, zmax
    double precision, parameter :: bohrToAngstrem = 0.529177249
    
    integer, dimension(:), allocatable :: nat, attract_data
    integer, dimension(:,:), allocatable :: index_basin
    integer*2, dimension(:), allocatable :: attract_code, attract_code_copy
    
    real, dimension(:), allocatable :: elf
    
    double precision, dimension(:), allocatable :: xat, yat, zat, x, y, z, volume
    
    character*4, dimension(:), allocatable :: atom_name
    character*16, dimension(:), allocatable :: basin_name, name_basin
    character*16 :: bname
    character*200 :: fileebas, filesij, fileelf, filenew
    
    logical if_file, read_elf_name, shutdown
    
    write(*,'(/,a,/)') 'TopologyParser 1.6 - Michal Michalski'
    write(*,'(a)') 'Podaj nazwe pliku ebas:'
    read(*,'(a200)') fileebas
    inquire(file = fileebas, exist = if_file)
    if (.not.if_file) then
        write(*,'("Plik nie istnieje: ",a200)') fileebas
        stop
    end if
    open(1, file = fileebas, status = 'unknown', form = 'unformatted')
    read(1) npx, npy, npz
    read(1) xmin, xmax, ymin, ymax, zmin, zmax
    read(1) natoms, attract_label
    
    npxyz = npx * npy * npz
    
    allocate(nat(natoms), atom_name(natoms), xat(natoms), yat(natoms), zat(natoms))
    do i = 1, natoms
        read(1) nat(i), atom_name(i), xat(i), yat(i), zat(i)
    end do
    
    allocate(volume(attract_label))
    allocate(index_basin(attract_label, 5), basin_name(attract_label), x(attract_label), y(attract_label), z(attract_label))
    do i = 1, attract_label
        read(1) index_basin(i, 1), basin_name(i), volume(i), x(i), y(i), z(i)
    end do
    
    allocate(attract_code(npxyz))
    read(1)(attract_code(i), i = 1, npxyz)
    close(1)
	
    write(*,'(a)') 'Podaj nazwe pliku esij:'
    read(*,'(a200)') filesij
	inquire(file = filesij, exist = if_file)
    if (if_file) then
        open(2, file = filesij, status='unknown', form='unformatted')
        read(2) nbas, nmol, natoms
        allocate(name_basin(nmol))
        read(2)(name_basin(i), i = 1, nbas)
        close(2)
        write(*,'(a,/)') 'TopologyParser bedzie wyswietlal nazwy basenow z populacji'
		read_elf_name = .true.
    else
        write(*,'(a,/)') 'TopologyParser bedzie wyswietlal nazwy basenow z listy atraktorow'
        read_elf_name = .false.
    end if
    
    write(*,'(a)') '1 - Zapisz atomy i atraktory do .xyz'
    write(*,'(a)') '2 - Modyfikuj plik ebas'
    write(*,'(a)') '3 - Modyfikuj plik elf'
    read(*,'(i1)') jobType
    
    allocate(attract_data(attract_label))
    do i = 1, size(attract_data)
        attract_data(i) = 0
    end do
    
    select case (jobType)
        case (1)
            filenew = fileebas(1 : (len_trim(fileebas) - 4))//'.xyz'
            open(2, file = filenew, status = 'unknown')
            write(2,*) natoms
            write(2,'(a)') 'TopologyParser'
            do i = 1, natoms
                write(2,'(i3,3(1x,f15.6))') nat(i), xat(i) * bohrToAngstrem, yat(i) * bohrToAngstrem, zat(i) * bohrToAngstrem
            end do
            close(2)
                
            filenew = fileebas(1 : (len_trim(fileebas) - 4))//'_att.xyz'
            open(3, file = filenew, status = 'unknown')
            write(3,*) attract_label
            write(3,'(a)') 'TopologyParser' 
            do i = 1, attract_label
                if (basin_name(i)(1:1).eq.'V') then
                    if (index(basin_name(i), ',') > 0) then
                        basin_name(i) = 'He'    ! V(N1,N2) or V(N1,N2,N3)
                    else
                        basin_name(i) = 'Ne'    ! V(N1)
                    end if
                else
                    basin_name(i) = 'Ar'        ! C(N1)
                end if
                write(3,'(a3,3(1x,f15.6))') basin_name(i), x(i) * bohrToAngstrem, y(i) * bohrToAngstrem, z(i) * bohrToAngstrem
            end do
            close(3)
        case (2)
            shutdown = .false.
            do while (.not.shutdown)
                write(*,'(a)') '1 - Usun baseny (grupowo)'
                write(*,'(a)') '2 - Usun basen (jeden)'
                write(*,'(a)') '3 - Usun wszystkie baseny poza wybranym'
                write(*,'(a)') '4 - Przyporzadkuj basenom grupe o nowym indeksie'
                write(*,'(a)') '5 - Wyswietl liste basenow'
                write(*,'(a)') '6 - Zapisz kazdy basen osobno'
                write(*,'(a)') '7 - Zapisz nowy plik ebas'
                read(*,'(i1)') jobType
				
                if (jobType == 4) then
                    write(*,'(a)') 'Podaj nowy indeks grupy'
                    read(*,'(i6)', iostat = iostat) newIndex
                end if
                
                if (jobType < 6) then
                    if (jobType == 1) then
                        do i = 1, attract_label
                            if (.not.read_elf_name) then
                                write(*,'(i6,1x,a,1x,a)') index_basin(i, 1), '---', basin_name(i)
                            else
                                write(*,'(i6,1x,a,1x,a)') index_basin(i, 1), '---', name_basin(i)
                            end if
                        end do
                    else
                        do i = 1, attract_label
                            if (.not.read_elf_name) then
                                write(*,'(i6,1x,a,1x,a)') indexList, '---', basin_name(i)
                            else
                                write(*,'(i6,1x,a,1x,a)') indexList, '---', name_basin(i)
                            end if
                            indexList = indexList + 1
                        end do
                    end if
                end if
            
                select case (jobType)
                    case (1, 2, 3, 4)
                        select case (jobType)
                            case (1)                           
                                write(*,'(a)') 'Wprowadz indeks basenu z listy:'
                                read(*,'(i6)', iostat = iostat) selection
                                if (iostat == 0) then
                                    do i = 1, attract_label
                                        if (index_basin(i, 1) == selection) then
                                            attract_data(i) = i
                                        end if
                                    end do
                                end if
                            case (2)
                                write(*,'(a)') 'Wprowadz numer basenu z listy:'
                                read(*,'(i6)', iostat = iostat) selection
                                if (iostat == 0) then
                                    attract_data(1) = selection
                                end if
                            case (3)                            
                                write(*,'(a)') 'Ile basenow zostawic?'
                                read(*,'(i6)', iostat = iostat) num
                                if (iostat == 0) then
                                    do i = 1, num
                                        write(*,'(a)') 'Wprowadz numer basenu z listy:'
                                        read(*,'(i6)', iostat = iostat) selection
                                        if (iostat == 0) then
                                            attract_data(i) = selection
                                        end if
                                    end do
                                end if
                            case (4)
								write(*,'(a)') 'Ile basenow zgrupowac?'
								read(*,'(i6)', iostat = iostat) num
								if (iostat == 0) then
									do i = 1, num
										write(*,'(a)') 'Wprowadz numer basenu z listy:'
										read(*,'(i6)', iostat = iostat) selection
										if (iostat == 0) then
											attract_data(i) = selection
										end if
									end do
								end if
                        end select

                        if (iostat == 0) then
                            do i = 1, npxyz
                                select case (jobType)
                                    case (1, 2)
                                        if (any(attract_data == attract_code(i))) then
                                            attract_code(i) = 0
                                        end if
                                    case (3)
                                        if (all(attract_data /= attract_code(i))) then
                                            attract_code(i) = 0
                                        end if
                                    case (4)
										if (any(attract_data == attract_code(i))) then
											attract_code(i) = newIndex
										end if
                                end select
                            end do
                        end if
                        
                        do i = 1, size(attract_data)
                            attract_data(i) = 0
                        end do
                        
                        indexList = 1
                        shutdown = .false.
                    case (5)
                        do i = 1, size(attract_data)
                            attract_data(i) = 0
                        end do
                        
                        indexList = 1
                        shutdown = .false.
                    case (6)
                        allocate(attract_code_copy(npxyz))
                        attract_code_copy = attract_code
                        
                        do i = 1, attract_label
                            attract_code = attract_code_copy
                            
                            do j = 1, npxyz
                                if (i /= attract_code(j)) then
                                    attract_code(j) = 0
                                end if
                            end do
                            
                            if (i < 10) then
                                write(filenew,'(a,i1,a)') '00', i, '.sbf'
                            else if (i < 100) then
                                write(filenew,'(a,i2,a)') '0', i, '.sbf'
                            else if (i < 1000) then
                                write(filenew,'(i3,a)') i, '.sbf'
                            end if

                            filenew = fileebas(1 : (len_trim(fileebas) - 4))//'_'//filenew
                            write(*,'(a,1x,a)') 'Zapisuje plik ---', filenew
                            
                            open(3, file = filenew, status = 'unknown', form = 'unformatted')
                            write(3) npx, npy, npz
                            write(3) xmin, xmax, ymin, ymax, zmin, zmax 
                            write(3) natoms, attract_label
                            do j = 1, natoms
                                write(3) nat(j), atom_name(j), xat(j), yat(j), zat(j)
                            end do
                            do j = 1, attract_label
                                write(3) index_basin(j, 1), basin_name(j), volume(j), x(j), y(j), z(j)
                            end do
                            write(3) attract_code
                            close(3)
                        end do
                        
                        shutdown = .true.
                    case (7)
                        filenew = fileebas(1 : (len_trim(fileebas) - 4))//'_NEW.sbf'
                        open(3, file = filenew, status = 'unknown', form = 'unformatted')
                        write(3) npx, npy, npz
                        write(3) xmin, xmax, ymin, ymax, zmin, zmax 
                        write(3) natoms, attract_label
                        do i = 1, natoms
                            write(3) nat(i), atom_name(i), xat(i), yat(i), zat(i)
                        end do
                        do i = 1, attract_label
                            write(3) index_basin(i, 1), basin_name(i), volume(i), x(i), y(i), z(i)
                        end do
                        write(3) attract_code
                        close(3)
                        shutdown = .true.
                end select
            end do          
        case (3)
            write(*,'(a)') 'Podaj nazwe pliku elf:'
            read(*,'(a200)') fileelf
            inquire(file = fileelf, exist = if_file)
            if(.not.if_file) then
                write(*,'("Plik nie istnieje: ",a200)') fileelf
                stop
            end if
            open(2, file = fileelf, status = 'unknown', form = 'unformatted')
            read(2) npx, npy, npz
            read(2) xmin, xmax, ymin, ymax, zmin, zmax
            
            npxyz = npx * npy * npz
            
            allocate(elf(npxyz))
            read(2)(elf(i), i = 1, npxyz)
            close(2)
            
            shutdown = .false.
            do while (.not.shutdown)
                write(*,'(a)') '1 - Usun domeny w obrebie basenu (grupowo)'
                write(*,'(a)') '2 - Usun domeny w obrebie basenu (wybrane)'
                write(*,'(a)') '3 - Usun wszystkie domeny w obrebie basenu poza wybranym'
                write(*,'(a)') '4 - Zapisz nowy plik ELF'
                read(*,'(i1)') jobType
				
                if (jobType < 4) then
                    if (jobType == 1) then
                        do i = 1, attract_label
                            if (.not.read_elf_name) then
                                write(*,'(i6,1x,a,1x,a)') index_basin(i, 1), '---', basin_name(i)
                            else
                                write(*,'(i6,1x,a,1x,a)') index_basin(i, 1), '---', name_basin(i)
                            end if
                        end do
                    else
                        do i = 1, attract_label
                            if (.not.read_elf_name) then
                                write(*,'(i6,1x,a,1x,a)') indexList, '---', basin_name(i)
                            else
                                write(*,'(i6,1x,a,1x,a)') indexList, '---', name_basin(i)
                            end if
                            indexList = indexList + 1
                        end do
                    end if
                end if
            
                select case (jobType)
                    case (1, 2, 3)
                        select case (jobType)
                            case (1)                           
                                write(*,'(a)') 'Wprowadz indeks basenu z listy:'
                                read(*,'(i6)', iostat = iostat) selection
                                if (iostat == 0) then
                                    do i = 1, attract_label
                                        if (index_basin(i, 1) == selection) then
                                            attract_data(i) = i
                                        end if
                                    end do
                                end if
                            case (2)
                                write(*,'(a)') 'Ile basenow usunac?'
                                read(*,'(i6)', iostat = iostat) num
                                if (iostat == 0) then
                                    do i = 1, num
                                        write(*,'(a)') 'Wprowadz numer basenu z listy:'
                                        read(*,'(i6)', iostat = iostat) selection
                                        if (iostat == 0) then
                                            attract_data(i) = selection
                                        end if
                                    end do
                                end if
                            case (3)
                                write(*,'(a)') 'Ile basenow zostawic?'
                                read(*,'(i6)', iostat = iostat) num
                                if (iostat == 0) then
                                    do i = 1, num
                                        write(*,'(a)') 'Wprowadz numer basenu z listy:'
                                        read(*,'(i6)', iostat = iostat) selection
                                        if (iostat == 0) then
                                            attract_data(i) = selection
                                        end if
                                    end do
                                end if
                        end select

                        if (iostat == 0) then
                            do i = 1, npxyz
                                select case (jobType)
                                    case (1, 2)
                                        if (any(attract_data == attract_code(i))) then
                                            elf(i) = 0
                                        end if
                                    case (3)
                                        if (all(attract_data /= attract_code(i))) then
                                            elf(i) = 0
                                        end if            
                                end select
                            end do
                        end if
                        
                        do i = 1, size(attract_data)
                            attract_data(i) = 0
                        end do
                        
                        indexList = 1
                        shutdown = .false.
                    case (4)
                        filenew = fileelf(1 : (len_trim(fileelf) - 4))//'_NEW.sbf'
                        open(3, file = filenew, status = 'unknown', form = 'unformatted')
                        write(3) npx, npy, npz
                        write(3) xmin, xmax, ymin, ymax, zmin, zmax 
                        write(3) elf
                        close(3)
                        shutdown = .true.
                end select
            end do
    end select
end program TopologyParser