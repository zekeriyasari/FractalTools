using SparseArray

struct Person
    name :: AbstractString
    id :: Int
end

people = Person[]

push!(people, Person("Gizem", 1001))