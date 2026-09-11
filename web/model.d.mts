export type State = [number, number, number, number];
export interface Parameters { m:number; I:number; delta:number; eps:number; k?:number; z0?:number; theta0?:number; }
export interface Energy { verticalEnergy:number; twistingEnergy:number; couplingEnergy:number; totalEnergy:number; }
export function validate(p:Parameters):{k:number;omega0Squared:number;splitting:number};
export function normalModes(p:Parameters):{plus:number;minus:number;natural:number};
export function analyticState(p:Parameters,elapsed:number,u0?:State):State;
export function derivative(p:Parameters,u:State):State;
export function rk4Step(p:Parameters,u:State,h:number):State;
export function energyOf(p:Parameters,u:State):Energy;
export function advanceClock(remainder:number,elapsed:number,speed:number,h:number,maxSteps?:number):{steps:number;remainder:number};
