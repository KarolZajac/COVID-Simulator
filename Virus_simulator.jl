#using Pkg
#Pkg.add("DynamicGrids")
#Pkg.add("ImageMagick")
#Pkg.instantiate()
#Pkg.add("Images")
#Pkg.add("FileIO")
#Pkg.add("ImageIO")

using DynamicGrids, ColorSchemes, Colors, Plots
using ImageView
using Images
using Pkg

SIZE = 400
ITERATIONS = 400
DEAD, RECOVERED, SUSCEPTIBLE, INFECTED = 1, 2, 3, 4
QUARANTINE_DAYS = 10
ANTIBODIES_PRESENCE_PERIOD = 90
#DEAD
#- RIP.
#RECOVERED:
#- State after INFECTION and spending QUARANTINE_DAYS,
#- Can be reinfected with very small reinfection_probability,
#- No symptoms, so can't die by virus
#- Resistant for virus for ANTIBODIES_PRESENCE_PERIOD and then become SUSCEPTIBLE

#SUSCEPTIBLE:
#- suspectible for infection, 
#- can't die by virus,
#- can get infection from neighbor with infect_probability,
#- can get infection by spawning infected cells with spawn_infected_probability

#INFECTED:
#- spending QUARANTINE_DAYS,
#- after QUARANTINE_DAYS become RECOVERED becouse cell get antibodies,
#- can infect SUSCEPTIBLE neighbors with infect_probability and RECOVERED neighbors with reinfection_probability,
#- can die with death_probability in health complications


A = zeros(ITERATIONS,4)
A[1,3]=SIZE*SIZE
quarantine_days = zeros(SIZE,SIZE)
antibodies_presence_days = zeros(SIZE,SIZE)

k=1
counter=0
new_cases= zeros(ITERATIONS)
total_cases = zeros(ITERATIONS)

spreadCOVID = let infect_probability=0.33, spawn_infected_probability=0.000001,
    death_probability=0.00001,reinfection_probability=0.000001
    SetNeighbors(Moore(1)) do data, neighborhood, cell, I
        if currentframe(data)==k+1
            global counter+=new_cases[k]
            global total_cases[k]=counter
            global k +=1
        end
        if cell == DEAD
            A[k,1]=A[k,1]+1 
        elseif cell == RECOVERED
            if antibodies_presence_days[I[1],I[2]] == ANTIBODIES_PRESENCE_PERIOD
                data[I...] = SUSCEPTIBLE
                antibodies_presence_days[I[1],I[2]]=0
            else
                antibodies_presence_days[I[1],I[2]]+=1
                A[k,2]=A[k,2]+1
            end 
        elseif cell == INFECTED
            if quarantine_days[I[1],I[2]] == QUARANTINE_DAYS
                data[I...] = RECOVERED
                quarantine_days[I[1],I[2]]=0
            else
                quarantine_days[I[1],I[2]]+=1
            end
            for pos in positions(neighborhood, I)
                if data[pos...] == SUSCEPTIBLE
                    if rand() <= infect_probability 
                        data[pos...] = INFECTED
                        new_cases[k]+=1
                    end
                elseif data[pos...] == RECOVERED
                    if rand() <= reinfection_probability 
                        data[pos...] = INFECTED
                        new_cases[k]+=1
                        antibodies_presence_days[pos[1],pos[2]]=0 #last edit
                    end
                end
            end
            if rand() <= death_probability
                data[I...] = DEAD
                A[k,1]=A[k,1]+1
                quarantine_days[I[1],I[2]]=0
                antibodies_presence_days[I[1],I[2]]=0
            else
                A[k,4]=A[k,4]+1 
            end
        elseif cell == SUSCEPTIBLE
            if rand() <= spawn_infected_probability 
                data[I...] = INFECTED
                A[k,4]=A[k,4]+1 
            else
                A[k,3]=A[k,3]+1
            end
        end
    end
end

init = fill(SUSCEPTIBLE, SIZE, SIZE)
output = GifOutput(init; 
    filename="./spread_virus.gif", tspan=1:ITERATIONS, fps=10, 
    minval=DEAD, maxval=INFECTED, 
    imagegen=Image(scheme=ColorSchemes.rainbow, zerocolor=RGB24(0.0)),
)

sim!(output, spreadCOVID)

state_change = @animate for i in 1:ITERATIONS
    plot(A[1:i,:],label=["DEAD" "RECOVERED" "SUSCEPRIBLE" "INFECTED"], xlabel="Time [Day]", ylabel="Population")
end

gif(state_change,"./state_change.gif",fps=10, variable_palette=true)

new_cases_plot = @animate for i in 1:ITERATIONS
    plot(new_cases[1:i],label="new_cases", xlabel="Time [Day]", ylabel="New cases")
end

gif(new_cases_plot,"./new_cases_plot.gif",fps=10, variable_palette=true)

total_cases_plot = @animate for i in 1:ITERATIONS
    plot(total_cases[1:i],label="total_cases", xlabel="Time [Day]", ylabel="Total cases")
end

gif(total_cases_plot,"./total_cases_plot.gif",fps=10, variable_palette=true)

using Images
using Pkg
#Pkg.add("ImageView")
using ImageView

visualisation = load("spread_virus.gif")
state_plot = load("state_change.gif")
imshow(visualisation)
imshow(state_plot)

